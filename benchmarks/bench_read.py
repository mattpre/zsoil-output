"""Benchmark (and cross-check) zsoil_output's binary result readers.

Run the same script against two checkouts of the package, then compare:

    git worktree add ../zsoil-output-main main
    python benchmarks/bench_read.py <project> --src ../zsoil-output-main --out main.npz
    python benchmarks/bench_read.py <project> --src .                     --out perf.npz
    python benchmarks/bench_read.py --compare main.npz perf.npz

<project> is the results path without extension (e.g. C:/runs/my-model/MY-MODEL).
"""
import argparse
import importlib
import os
import statistics
import subprocess
import sys
import time

import numpy as np


def import_package(src):
    src = os.path.abspath(src)
    if os.path.isdir(os.path.join(src, "src", "zsoil_output")):
        src = os.path.join(src, "src")
    if not os.path.isdir(os.path.join(src, "zsoil_output")):
        sys.exit("no zsoil_output package found under %s" % src)
    sys.path.insert(0, src)
    pkg = importlib.import_module("zsoil_output")
    mods = {}
    for name in ["C_Mesh", "C_HistoryOfExecution", "C_Rcf_info", "C_Rcf_layers_info",
                 "C_NodalResults", "C_EleResults", "C_EleLayersResults"]:
        mods[name] = importlib.import_module("zsoil_output." + name)
    try:
        rev = subprocess.run(["git", "-C", src, "log", "-1", "--format=%h %s"],
                             capture_output=True, text=True).stdout.strip()
    except OSError:
        rev = "?"
    print("zsoil_output from %s  [%s]" % (os.path.dirname(pkg.__file__), rev))
    return mods


def subsample(items, n_max):
    if n_max is None or len(items) <= n_max:
        return list(items)
    idx = np.linspace(0, len(items) - 1, n_max).round().astype(int)
    return [items[i] for i in np.unique(idx)]


def timed(fn, repeat):
    times = []
    result = None
    for _ in range(repeat):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return result, times


def as_array(res):
    if isinstance(res, tuple):
        return np.stack([np.asarray(r, dtype=np.float64) for r in res])
    if isinstance(res, list):
        if len(res) == 0:
            return np.zeros(0)
        return np.asarray([np.asarray(r, dtype=np.float64) for r in res])
    return np.asarray(res, dtype=np.float64)


def run(args):
    m = import_package(args.src)
    project = args.project

    t0 = time.perf_counter()
    his = m["C_HistoryOfExecution"].HistoryOfExecution(project)
    rcf = m["C_Rcf_info"].RCF_info(project)
    rcf_lay = m["C_Rcf_layers_info"].RCF_layers_info(project)
    mesh = m["C_Mesh"].Mesh(project)
    t_load = time.perf_counter() - t0

    sols = subsample(his.give_converged_solutions(), args.max_steps)
    node_ids = subsample([n.index for n in mesh.nodes], args.max_nodes)
    print("loaded mesh/his/rcf in %.2f s: %d nodes, %d elements, %d of %d converged steps used"
          % (t_load, len(mesh.nodes), len(mesh.elements), len(sols), len(his.give_converged_solutions())))

    work = []  # (name, callable)

    if os.path.exists(project + ".s00"):
        def nodal_all_steps():
            nd = m["C_NodalResults"].NodalResults(mesh, his, rcf)
            return [nd.get_nodes_displacements(node_ids, s) for s in sols]

        def nodal_time_history():
            nd = m["C_NodalResults"].NodalResults(mesh, his, rcf)
            return [nd.get_node_displacements_time_history(i, sols) for i in node_ids[:args.n_th_nodes]]

        def nodal_envelope():
            nd = m["C_NodalResults"].NodalResults(mesh, his, rcf)
            return nd.get_displ_envelope_for_sel_nodes(node_ids, sols, 'ABS')

        work += [("nodal_all_steps", nodal_all_steps),
                 ("nodal_time_history", nodal_time_history),
                 ("nodal_envelope", nodal_envelope)]

    for group in mesh.groups_with_rsl:
        elems = mesh.get_list_of_elements(group)
        if not elems or not os.path.exists(project + elems[0].rsl_ext):
            continue
        elems = subsample(elems, args.max_elements)

        def ele_raw(elems=elems):
            er = m["C_EleResults"].EleResults(mesh, his, rcf)
            return er.get_element_results_ex(elems, sols)

        def ele_one_by_one(elems=elems):
            # typical post-processing pattern: one element (all steps) per call
            er = m["C_EleResults"].EleResults(mesh, his, rcf)
            out = []
            for e in elems[:args.n_th_elements]:
                out += er.get_element_results_ex([e], sols)
            return out

        work += [("ele_%s" % group, ele_raw), ("ele_%s_one_by_one" % group, ele_one_by_one)]

    for group in ["BEAMS", "SHELLS"]:
        elems = mesh.get_list_of_elements(group)
        if not elems or elems[0].rsl_lay_ext is None or not os.path.exists(project + elems[0].rsl_lay_ext):
            continue

        def layers(group=group):
            lr = m["C_EleLayersResults"].EleLayersResults(mesh, his, rcf, rcf_lay)
            return lr.get_layers_results(group, sols)

        work.append(("layers_%s" % group, layers))

    results = {}
    print("%-28s %10s %10s   (seconds, %d repeats)" % ("workload", "min", "median", args.repeat))
    for name, fn in work:
        try:
            res, times = timed(fn, args.repeat)
        except Exception as exc:  # keep going so one broken workload doesn't hide the others
            print("%-28s FAILED: %r" % (name, exc))
            continue
        print("%-28s %10.4f %10.4f" % (name, min(times), statistics.median(times)))
        results[name] = as_array(res)
        results["time__" + name] = np.asarray(times)

    if args.out:
        np.savez_compressed(args.out, **results)
        print("saved results to %s" % args.out)


def compare(a_path, b_path):
    a = np.load(a_path)
    b = np.load(b_path)
    names = sorted(k for k in set(a.files) | set(b.files) if not k.startswith("time__"))
    ok = True
    print("%-28s %10s %10s %8s  %s" % ("workload", "A min[s]", "B min[s]", "speedup", "results"))
    for k in names:
        if k not in a.files or k not in b.files:
            print("%-28s missing in %s" % (k, "A" if k not in a.files else "B"))
            ok = False
            continue
        same = a[k].shape == b[k].shape and np.array_equal(a[k], b[k], equal_nan=True)
        ok = ok and same
        ta, tb = a["time__" + k].min(), b["time__" + k].min()
        print("%-28s %10.4f %10.4f %7.1fx  %s" % (k, ta, tb, ta / tb if tb > 0 else float("inf"),
                                                "identical" if same else "DIFFERENT"))
    print("all results identical" if ok else "RESULTS DIFFER")
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("project", nargs="?", help="results path without extension")
    p.add_argument("--src", default=os.path.join(os.path.dirname(__file__), ".."),
                   help="checkout (or its src/ dir) to import zsoil_output from")
    p.add_argument("--out", help="save results + timings to this .npz")
    p.add_argument("--compare", nargs=2, metavar=("A.npz", "B.npz"))
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--max-steps", type=int, default=None, help="subsample converged steps")
    p.add_argument("--max-nodes", type=int, default=None, help="subsample nodes")
    p.add_argument("--max-elements", type=int, default=5000, help="subsample elements per group")
    p.add_argument("--n-th-nodes", type=int, default=50, help="nodes for the time-history workload")
    p.add_argument("--n-th-elements", type=int, default=200, help="elements for the one-by-one workload")
    args = p.parse_args()

    if args.compare:
        sys.exit(compare(*args.compare))
    if not args.project:
        p.error("project is required unless --compare is given")
    run(args)


if __name__ == "__main__":
    main()
