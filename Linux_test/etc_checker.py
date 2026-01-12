import difflib


def normalize_output(output):
    if output is None:
        return []
    lines = output.splitlines()
    normalized = []
    for line in lines:
        stripped = line.strip()
        if stripped:
            normalized.append(stripped)
    return normalized


def compare_outputs(rtl_output, netlist_output):
    rtl_lines = normalize_output(rtl_output)
    netlist_lines = normalize_output(netlist_output)
    if rtl_lines == netlist_lines:
        return True, ""

    diff = difflib.unified_diff(
        rtl_lines,
        netlist_lines,
        fromfile="rtl",
        tofile="netlist",
        lineterm="",
    )
    return False, "\n".join(diff)
