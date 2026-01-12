classdef StateflowTPG < emi.decs.DecoratedMutator
    %STATEFLOWTPG Apply Stateflow-based structural mutations (TPG).

    methods
        function obj = StateflowTPG(varargin)
            obj = obj@emi.decs.DecoratedMutator(varargin{:});
        end

        function main_phase(obj)
            if ~emi.cfg.TPG_ENABLE
                return;
            end

            charts = obj.find_stateflow_charts();
            if isempty(charts)
                obj.l.info('TPG: No Stateflow charts found in %s', obj.mutant.sys);
                return;
            end

            for i = 1:numel(charts)
                chart = charts{i};
                obj.ensure_selector_data(chart);

                for j = 1:emi.cfg.TPG_MUTATIONS_PER_CHART
                    strategy_id = randi(3);
                    applied = false;

                    switch strategy_id
                        case 1
                            applied = obj.state_duplication(chart);
                        case 2
                            applied = obj.path_duplication(chart);
                        case 3
                            applied = obj.transition_expansion(chart);
                    end

                    if applied
                        obj.l.info('TPG: Applied strategy %d on chart %s', strategy_id, chart.Path);
                    end
                end
            end
        end
    end

    methods(Access = private)
        function charts = find_stateflow_charts(obj)
            rt = sfroot;
            all_charts = rt.find('-isa', 'Stateflow.Chart');
            charts = {};
            if isempty(all_charts)
                return;
            end

            prefix = [obj.mutant.sys '/'];
            for i = 1:numel(all_charts)
                chart = all_charts(i);
                if strcmp(chart.Path, obj.mutant.sys) || strncmp(chart.Path, prefix, length(prefix))
                    charts{end + 1} = chart; %#ok<AGROW>
                end
            end
        end

        function ensure_selector_data(~, chart)
            selector_name = emi.cfg.TPG_SELECTOR_NAME;
            existing = chart.find('-isa', 'Stateflow.Data', 'Name', selector_name);
            if ~isempty(existing)
                return;
            end

            selector = Stateflow.Data(chart);
            selector.Name = selector_name;
            selector.Scope = 'Local';
            selector.DataType = 'boolean';
            if isprop(selector, 'Props')
                selector.Props.InitialValue = '0';
            end
        end

        function applied = state_duplication(obj, chart)
            applied = false;
            states = chart.find('-isa', 'Stateflow.State');
            if isempty(states)
                return;
            end

            state = states(randi(numel(states)));
            if obj.has_substates(state)
                return;
            end

            clone = obj.clone_state(chart, state, '_dup');
            if isempty(clone)
                return;
            end

            transitions = chart.find('-isa', 'Stateflow.Transition');
            if isempty(transitions)
                return;
            end

            for i = 1:numel(transitions)
                t = transitions(i);
                if isempty(t.Destination) || isempty(t.Source)
                    continue;
                end

                if t.Destination == state
                    original_label = t.LabelString;
                    t.LabelString = obj.add_guard_to_label(original_label, ['~' emi.cfg.TPG_SELECTOR_NAME]);

                    new_t = Stateflow.Transition(chart);
                    new_t.Source = t.Source;
                    new_t.Destination = clone;
                    new_t.LabelString = obj.add_guard_to_label(original_label, emi.cfg.TPG_SELECTOR_NAME);
                elseif t.Source == state
                    new_t = Stateflow.Transition(chart);
                    new_t.Source = clone;
                    new_t.Destination = t.Destination;
                    new_t.LabelString = t.LabelString;
                end
            end

            applied = true;
        end

        function applied = path_duplication(obj, chart)
            applied = false;
            states = chart.find('-isa', 'Stateflow.State');
            transitions = chart.find('-isa', 'Stateflow.Transition');
            if isempty(states) || isempty(transitions)
                return;
            end

            [indeg, outdeg] = obj.compute_degrees(states, transitions);
            eligible = states(indeg == 1 & outdeg == 1);
            if isempty(eligible)
                return;
            end

            mid = eligible(randi(numel(eligible)));
            if obj.has_substates(mid)
                return;
            end

            chain = obj.expand_linear_chain(mid, states, transitions, indeg, outdeg);
            if numel(chain) < 2
                return;
            end

            chain_ids = arrayfun(@(s) s.Id, chain);
            entry_transitions = obj.transitions_into_state(transitions, chain(1));
            entry_transitions = entry_transitions(arrayfun(@(t) ~obj.is_state_in_chain(t.Source, chain_ids), entry_transitions));
            if isempty(entry_transitions)
                return;
            end
            entry_t = entry_transitions(1);

            exit_transitions = obj.transitions_out_of_state(transitions, chain(end));
            exit_transitions = exit_transitions(arrayfun(@(t) ~obj.is_state_in_chain(t.Destination, chain_ids), exit_transitions));
            if isempty(exit_transitions)
                return;
            end
            exit_t = exit_transitions(1);

            clone_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:numel(chain)
                clone_map(chain(i).Id) = obj.clone_state(chart, chain(i), '_dup');
            end

            original_label = entry_t.LabelString;
            entry_t.LabelString = obj.add_guard_to_label(original_label, ['~' emi.cfg.TPG_SELECTOR_NAME]);
            new_entry = Stateflow.Transition(chart);
            new_entry.Source = entry_t.Source;
            new_entry.Destination = clone_map(chain(1).Id);
            new_entry.LabelString = obj.add_guard_to_label(original_label, emi.cfg.TPG_SELECTOR_NAME);

            for i = 1:numel(transitions)
                t = transitions(i);
                if isempty(t.Source) || isempty(t.Destination)
                    continue;
                end
                if obj.is_state_in_chain(t.Source, chain_ids) && obj.is_state_in_chain(t.Destination, chain_ids)
                    new_t = Stateflow.Transition(chart);
                    new_t.Source = clone_map(t.Source.Id);
                    new_t.Destination = clone_map(t.Destination.Id);
                    new_t.LabelString = t.LabelString;
                end
            end

            new_exit = Stateflow.Transition(chart);
            new_exit.Source = clone_map(chain(end).Id);
            new_exit.Destination = exit_t.Destination;
            new_exit.LabelString = exit_t.LabelString;

            applied = true;
        end

        function applied = transition_expansion(~, chart)
            applied = false;
            transitions = chart.find('-isa', 'Stateflow.Transition');
            if isempty(transitions)
                return;
            end

            t = transitions(randi(numel(transitions)));
            if isempty(t.Source) || isempty(t.Destination)
                return;
            end

            old_dest = t.Destination;
            junction = Stateflow.Junction(chart);
            junction.Position = StateflowTPG.midpoint(t.Source, old_dest);

            t.Destination = junction;

            new_t = Stateflow.Transition(chart);
            new_t.Source = junction;
            new_t.Destination = old_dest;
            new_t.LabelString = '[true]';

            applied = true;
        end

        function [indeg, outdeg] = compute_degrees(~, states, transitions)
            indeg = zeros(size(states));
            outdeg = zeros(size(states));
            state_ids = arrayfun(@(s) s.Id, states);
            for i = 1:numel(transitions)
                t = transitions(i);
                if isempty(t.Source) || isempty(t.Destination)
                    continue;
                end
                if isa(t.Source, 'Stateflow.State') && isa(t.Destination, 'Stateflow.State')
                    src_idx = find(state_ids == t.Source.Id, 1);
                    dst_idx = find(state_ids == t.Destination.Id, 1);
                    if ~isempty(src_idx)
                        outdeg(src_idx) = outdeg(src_idx) + 1;
                    end
                    if ~isempty(dst_idx)
                        indeg(dst_idx) = indeg(dst_idx) + 1;
                    end
                end
            end
        end

        function chain = expand_linear_chain(obj, start_state, states, transitions, indeg, outdeg)
            chain = start_state;

            state_ids = arrayfun(@(s) s.Id, states);
            start_idx = find(state_ids == start_state.Id, 1);
            if isempty(start_idx)
                return;
            end

            current = start_state;
            while true
                next_transitions = obj.transitions_out_of_state(transitions, current);
                if numel(next_transitions) ~= 1
                    break;
                end
                next_state = next_transitions(1).Destination;
                if isempty(next_state) || ~isa(next_state, 'Stateflow.State')
                    break;
                end
                next_idx = find(state_ids == next_state.Id, 1);
                if isempty(next_idx) || indeg(next_idx) ~= 1 || outdeg(next_idx) ~= 1
                    break;
                end
                if obj.has_substates(next_state)
                    break;
                end
                chain(end + 1) = next_state; %#ok<AGROW>
                current = next_state;
            end
        end

        function transitions = transitions_into_state(~, transitions, state)
            transitions = transitions(arrayfun(@(t) ~isempty(t.Destination) && t.Destination == state, transitions));
        end

        function transitions = transitions_out_of_state(~, transitions, state)
            transitions = transitions(arrayfun(@(t) ~isempty(t.Source) && t.Source == state, transitions));
        end

        function ret = is_state_in_chain(~, obj_state, chain_ids)
            ret = isa(obj_state, 'Stateflow.State') && any(chain_ids == obj_state.Id);
        end

        function ret = has_substates(~, state)
            ret = false;
            if isprop(state, 'Substates')
                ret = ~isempty(state.Substates);
            end
        end

        function clone = clone_state(~, chart, state, suffix)
            clone = Stateflow.State(chart);
            clone.Name = [state.Name suffix];
            if isprop(state, 'LabelString')
                clone.LabelString = state.LabelString;
            end
            if isprop(state, 'Position')
                pos = state.Position;
                if numel(pos) >= 2
                    pos(1) = pos(1) + 80;
                    pos(2) = pos(2) + 80;
                end
                clone.Position = pos;
            end
        end

        function new_label = add_guard_to_label(obj, label, guard)
            [event, existing_guard, action] = obj.parse_transition_label(label);
            if isempty(existing_guard)
                merged_guard = guard;
            else
                merged_guard = ['(' existing_guard ') && (' guard ')'];
            end
            new_label = obj.build_transition_label(event, merged_guard, action);
        end

        function [event, guard, action] = parse_transition_label(~, label)
            guard_match = regexp(label, '\[(.*?)\]', 'tokens', 'once');
            action_match = regexp(label, '\{(.*?)\}', 'tokens', 'once');

            guard = '';
            action = '';
            if ~isempty(guard_match)
                guard = strtrim(guard_match{1});
            end
            if ~isempty(action_match)
                action = strtrim(action_match{1});
            end

            event = strtrim(regexprep(label, '\[.*?\]|\{.*?\}', ''));
        end

        function label = build_transition_label(~, event, guard, action)
            label = strtrim(event);
            if ~isempty(guard)
                if isempty(label)
                    label = ['[' guard ']'];
                else
                    label = [label ' [' guard ']'];
                end
            end
            if ~isempty(action)
                if isempty(label)
                    label = ['{' action '}'];
                else
                    label = [label ' {' action '}'];
                end
            end
        end
    end

    methods(Static, Access = private)
        function pos = midpoint(src, dst)
            src_pos = StateflowTPG.object_position(src);
            dst_pos = StateflowTPG.object_position(dst);
            pos = (src_pos + dst_pos) ./ 2;
        end

        function pos = object_position(obj)
            pos = [0 0];
            if isprop(obj, 'Position')
                raw = obj.Position;
                if numel(raw) == 4
                    pos = [raw(1) + raw(3) / 2, raw(2) + raw(4) / 2];
                elseif numel(raw) >= 2
                    pos = raw(1:2);
                end
            end
        end
    end
end
