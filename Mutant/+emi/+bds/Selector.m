classdef Selector < handle
    %SELECTOR Bayesian Diversity Selection (BDS) for mutant candidates.

    methods
        function [selected_idx, scored_reports] = select(obj, reports)
            valid_idx = find(~cellfun(@(r) isempty(r) || ~isempty(r.exception), reports));
            if isempty(valid_idx)
                selected_idx = [];
                scored_reports = reports;
                return;
            end

            [features, pressures, scored_reports] = obj.score_reports(reports, valid_idx);
            scores = obj.combine_scores(features, pressures);

            if emi.cfg.BDS_USE_SMBO
                selected_idx = obj.smbo_select(valid_idx, features, scores);
            else
                selected_idx = obj.greedy_select(valid_idx, features, scores);
            end
        end
    end

    methods(Access = private)
        function [features, pressures, reports] = score_reports(obj, reports, valid_idx)
            num_valid = numel(valid_idx);
            features = zeros(num_valid, 6);
            pressures = zeros(num_valid, 3);

            for i = 1:num_valid
                idx = valid_idx(i);
                report = reports{idx};
                [feat_vec, graph_metrics] = obj.extract_stateflow_features(report);
                [pressure_vec, pressure_metrics] = obj.extract_pressure(report, graph_metrics);

                features(i, :) = feat_vec;
                pressures(i, :) = pressure_vec;

                report.bds_features = feat_vec;
                report.bds_pressure = pressure_vec;
                report.bds_graph = graph_metrics;
                report.bds_pressure_meta = pressure_metrics;
                reports{idx} = report;
            end
        end

        function scores = combine_scores(obj, features, pressures)
            norm_features = obj.minmax_normalize(features);
            norm_pressure = obj.minmax_normalize(pressures);

            pressure_score = min(norm_pressure, [], 2);

            diversity = obj.diversity_distance(norm_features);
            diversity = obj.minmax_normalize(diversity);

            scores = pressure_score .* diversity;
        end

        function selected_idx = greedy_select(obj, valid_idx, features, scores)
            num_select = min(emi.cfg.BDS_SELECT_COUNT, numel(valid_idx));
            selected_idx = zeros(1, num_select);
            norm_features = obj.minmax_normalize(features);

            remaining = valid_idx;
            for i = 1:num_select
                remaining_scores = scores(ismember(valid_idx, remaining));
                remaining_features = norm_features(ismember(valid_idx, remaining), :);

                if i == 1
                    [~, max_idx] = max(remaining_scores);
                    selected_idx(i) = remaining(max_idx);
                else
                    selected_features = norm_features(ismember(valid_idx, selected_idx(1:i-1)), :);
                    dists = obj.min_distance(remaining_features, selected_features);
                    [~, max_idx] = max(remaining_scores .* dists);
                    selected_idx(i) = remaining(max_idx);
                end

                remaining(remaining == selected_idx(i)) = [];
            end
        end

        function selected_idx = smbo_select(obj, valid_idx, features, scores)
            num_select = min(emi.cfg.BDS_SELECT_COUNT, numel(valid_idx));
            selected_idx = zeros(1, num_select);

            if ~exist('TreeBagger', 'class')
                selected_idx = obj.greedy_select(valid_idx, features, scores);
                return;
            end

            norm_features = obj.minmax_normalize(features);
            candidate_idx = 1:numel(valid_idx);

            [~, best_idx] = max(scores);
            selected_local = best_idx;

            for i = 2:num_select
                train_x = norm_features(selected_local, :);
                train_y = scores(selected_local);

                model = TreeBagger(...
                    emi.cfg.BDS_RF_TREES, train_x, train_y, ...
                    'Method', 'regression', ...
                    'MinLeafSize', emi.cfg.BDS_RF_MIN_LEAF ...
                );

                remaining_local = setdiff(candidate_idx, selected_local);
                test_x = norm_features(remaining_local, :);
                [pred_mean, pred_std] = obj.rf_predict(model, test_x);

                best_score = max(train_y);
                ei = obj.expected_improvement(pred_mean, pred_std, best_score);
                score = max(ei, pred_mean);
                [~, pick_idx] = max(score);
                selected_local(end + 1) = remaining_local(pick_idx); %#ok<AGROW>
            end

            selected_idx = valid_idx(selected_local);
        end

        function [features, metrics] = extract_stateflow_features(~, report)
            metrics = struct('nodes', 0, 'edges', 0, 'components', 0);
            features = zeros(1, 6);

            [loaded, mdl, loc] = emi.bds.Selector.load_model(report);
            if ~loaded
                emi.bds.Selector.cleanup_path(loc);
                return;
            end

            root = sfroot;
            charts = root.find('-isa', 'Stateflow.Chart', 'Path', mdl);
            if isempty(charts)
                emi.bds.Selector.close_model(mdl);
                emi.bds.Selector.cleanup_path(loc);
                return;
            end

            states = charts.find('-isa', 'Stateflow.State');
            transitions = charts.find('-isa', 'Stateflow.Transition');
            junctions = charts.find('-isa', 'Stateflow.Junction');

            nodes = [states; junctions];
            node_ids = arrayfun(@(n) n.Id, nodes);
            edges = zeros(numel(transitions), 2);
            edge_count = 0;

            outdeg = zeros(numel(nodes), 1);
            indeg = zeros(numel(nodes), 1);

            for i = 1:numel(transitions)
                t = transitions(i);
                if isempty(t.Source) || isempty(t.Destination)
                    continue;
                end
                if ~isa(t.Source, 'Stateflow.State') && ~isa(t.Source, 'Stateflow.Junction')
                    continue;
                end
                if ~isa(t.Destination, 'Stateflow.State') && ~isa(t.Destination, 'Stateflow.Junction')
                    continue;
                end
                src_idx = find(node_ids == t.Source.Id, 1);
                dst_idx = find(node_ids == t.Destination.Id, 1);
                if isempty(src_idx) || isempty(dst_idx)
                    continue;
                end
                edge_count = edge_count + 1;
                edges(edge_count, :) = [src_idx, dst_idx];
                outdeg(src_idx) = outdeg(src_idx) + 1;
                indeg(dst_idx) = indeg(dst_idx) + 1;
            end

            edges = edges(1:edge_count, :);

            diameter = emi.bds.Selector.graph_diameter(numel(nodes), edges);
            branch_nodes = sum(outdeg > 1);
            merge_nodes = sum(indeg > 1);

            features = [numel(states), numel(transitions), numel(junctions), diameter, branch_nodes, merge_nodes];

            metrics.nodes = numel(nodes);
            metrics.edges = edge_count;
            metrics.components = emi.bds.Selector.count_components(numel(nodes), edges);

            emi.bds.Selector.close_model(mdl);
            emi.bds.Selector.cleanup_path(loc);
        end

        function [pressure, metrics] = extract_pressure(obj, report, graph_metrics)
            metrics = struct('hdl_files', 0, 'bus_width', 0, 'fsm_complexity', 0, 'mux_density', 0);
            [bus_width, mux_density] = obj.scan_hdl_metrics(report);

            fsm_complexity = 0;
            if graph_metrics.nodes > 0
                fsm_complexity = max(0, graph_metrics.edges - graph_metrics.nodes + graph_metrics.components);
            end

            pressure = [bus_width, fsm_complexity, mux_density];
            metrics.bus_width = bus_width;
            metrics.fsm_complexity = fsm_complexity;
            metrics.mux_density = mux_density;
        end

        function [bus_width, mux_density] = scan_hdl_metrics(~, report)
            bus_width = 0;
            mux_density = 0;

            roots = emi.cfg.BDS_HDL_SEARCH_DIRS;
            if isempty(roots)
                roots = {report.loc};
            end

            all_files = [];
            for i = 1:numel(roots)
                if isempty(roots{i}) || ~isfolder(roots{i})
                    continue;
                end
                files = dir(fullfile(roots{i}, '**', '*.v'));
                all_files = [all_files; files]; %#ok<AGROW>
            end

            if isempty(all_files)
                return;
            end

            total_lines = 0;
            mux_count = 0;

            for i = 1:numel(all_files)
                path = fullfile(all_files(i).folder, all_files(i).name);
                content = fileread(path);

                widths = regexp(content, '\[(\d+)\s*:\s*(\d+)\]', 'tokens');
                for j = 1:numel(widths)
                    hi = str2double(widths{j}{1});
                    lo = str2double(widths{j}{2});
                    bus_width = bus_width + abs(hi - lo) + 1;
                end

                mux_count = mux_count + numel(regexp(content, '\?', 'match'));
                mux_count = mux_count + numel(regexp(content, '\bcase\b', 'match'));

                total_lines = total_lines + max(1, numel(strfind(content, newline)));
            end

            mux_density = mux_count / max(1, total_lines);
        end

        function normed = minmax_normalize(~, values)
            if isempty(values)
                normed = values;
                return;
            end
            mins = min(values, [], 1);
            maxs = max(values, [], 1);
            ranges = maxs - mins;
            ranges(ranges == 0) = 1;
            normed = (values - mins) ./ ranges;
        end

        function diversity = diversity_distance(obj, features)
            if isempty(features)
                diversity = [];
                return;
            end
            diversity = zeros(size(features, 1), 1);
            for i = 1:size(features, 1)
                if i == 1
                    diversity(i) = 1;
                    continue;
                end
                other = features(1:i-1, :);
                distances = obj.min_distance(features(i, :), other);
                diversity(i) = distances;
            end
        end

        function d = min_distance(~, points, refs)
            if isempty(refs)
                d = ones(size(points, 1), 1);
                return;
            end
            if size(points, 1) == 1
                distances = sqrt(sum((refs - points).^2, 2));
                d = min(distances);
            else
                d = zeros(size(points, 1), 1);
                for i = 1:size(points, 1)
                    distances = sqrt(sum((refs - points(i, :)).^2, 2));
                    d(i) = min(distances);
                end
            end
        end

        function [mu, sigma] = rf_predict(~, model, x)
            preds = predict(model, x);
            if iscell(preds)
                preds = cellfun(@str2double, preds);
            end

            if isvector(preds)
                mu = preds;
                sigma = zeros(size(mu));
                return;
            end

            mu = mean(preds, 2);
            sigma = std(preds, 0, 2);
        end

        function ei = expected_improvement(~, mu, sigma, best)
            sigma = max(sigma, 1e-6);
            z = (mu - best) ./ sigma;
            ei = (mu - best) .* normcdf(z) + sigma .* normpdf(z);
            ei = max(ei, 0);
        end
    end

    methods(Static, Access = private)
        function [loaded, model_name, loc] = load_model(report)
            loaded = false;
            model_name = report.sys;
            loc = report.loc;

            if isempty(report.loc)
                return;
            end

            addpath(report.loc);
            try
                load_system(model_name);
                loaded = true;
            catch
                loaded = false;
            end
        end

        function cleanup_path(loc)
            if ~isempty(loc)
                rmpath(loc);
            end
        end

        function close_model(model_name)
            if bdIsLoaded(model_name)
                bdclose(model_name);
            end
        end

        function diameter = graph_diameter(node_count, edges)
            if node_count == 0
                diameter = 0;
                return;
            end
            adjacency = cell(node_count, 1);
            for i = 1:size(edges, 1)
                src = edges(i, 1);
                dst = edges(i, 2);
                adjacency{src} = [adjacency{src}, dst];
                adjacency{dst} = [adjacency{dst}, src];
            end

            diameter = 0;
            for i = 1:node_count
                dist = emi.bds.Selector.bfs_distances(adjacency, i);
                diameter = max(diameter, max(dist(dist < inf)));
            end
        end

        function components = count_components(node_count, edges)
            if node_count == 0
                components = 0;
                return;
            end
            adjacency = cell(node_count, 1);
            for i = 1:size(edges, 1)
                src = edges(i, 1);
                dst = edges(i, 2);
                adjacency{src} = [adjacency{src}, dst];
                adjacency{dst} = [adjacency{dst}, src];
            end

            visited = false(node_count, 1);
            components = 0;
            for i = 1:node_count
                if visited(i)
                    continue;
                end
                components = components + 1;
                queue = i;
                visited(i) = true;
                while ~isempty(queue)
                    current = queue(1);
                    queue(1) = [];
                    neighbors = adjacency{current};
                    for j = 1:numel(neighbors)
                        n = neighbors(j);
                        if ~visited(n)
                            visited(n) = true;
                            queue(end + 1) = n; %#ok<AGROW>
                        end
                    end
                end
            end
        end

        function dist = bfs_distances(adjacency, start_node)
            node_count = numel(adjacency);
            dist = inf(node_count, 1);
            dist(start_node) = 0;
            queue = start_node;
            while ~isempty(queue)
                current = queue(1);
                queue(1) = [];
                neighbors = adjacency{current};
                for i = 1:numel(neighbors)
                    n = neighbors(i);
                    if dist(n) == inf
                        dist(n) = dist(current) + 1;
                        queue(end + 1) = n; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
