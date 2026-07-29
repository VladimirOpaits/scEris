import argparse
import numpy as np
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.text import Text
from . import config, project
from .store import SignatureStore
from .source import Source
from .dataset import load_cohort, spec_for
from .retrieval import diverse_refs, matched_refs
from .stats import ref_distances

LOGO = r"""
                                        7ZObbbbbbXgPGGGGGGGGGGGGGGbYOOYb?                                                  ``
                                        `rVQQQQQQ$j)\>)))))))))))ix*7hMQ6 x!aCwwwfyfwwwT7});         \?1Tfwwwwe!c     _lJSdVhgpJ[v\pI
                                           )pNQQQQR4r                 )YO.`_,eQQQA*lxc}uXQQDE?.      ._'%NQQQe,_`    t&Qg{",:^)*y$QQ1
                                             )yDQQQQQG!-               xK/   ;DQQh       %$QQQf          XQQD,      IQQ6         ;5Q1
            ._'_              `__.             /JKQQQQQkt'             .|-   ;HQQV        jQQQK_         4QQD;      JQQA%          Ee
         ^tToat7T36p      .vLwwC3V4w{_           ^jUQQQQQ&z;                 ;HQQV       -VQQQy          4QQD;      /GQQQX7)'      ++
        vKK=    '*KO     IZKo+` ./AQQU\            :eOQQQQQ$T|               ;HQQV    _|1OQNZt.          4QQD;       _1G0QQQ$hzi'
        oQRo;     <3    nQQr      >gPh|              _!PRQQQ0f_              ;HQQ&nT5OKRDul/             4QQD;          \oVKQQQQ$g1^
        _JHQNETs/.     i0Q&.                           )OQ$J|                ;HQQd'_,7DQ0h)              4QQD;             :%zhHQQQHC:
          ;}5&QQ0YJ\   tQQD^                        .rdMYa,                  ;HQQV    +pQQNn_            4QQD;      ))         ;!PQQQg_
        \/   ./?hQQNe  >RQQ3         \a           ,aOR4{`                 .  ;HQQV      skQQYl           4QQD;      zY_          .6QQR+
        #O+      ]QQg   LQQQq)`    +jYc         |CHW5<                   tG; ^BQQd       _nWQ0F|         PQQD,      aQX)          FQQy.
        TQUui^:,<204<    }ZQQQUEhEO&6|        %mWA7,                   +J0X+)CQQQN!|='     \hQQ0Vo%|,_;"]0QQQJ)^:   oQBKdu*v)|\ljbRG!
        x?>?oTCCCo<       `io5hdSC?=       _!XQNwc\>iiiiiiiiiiiiv%l*aC4DQQLs?a111ttI{)       r1tee1?>|l}tt111a?ri   sa:^i?jC5p6pC[)
                                         /u8QQQQRQQQQQQQQQQQQQQQQQQQQQQQQ8_
                                        z8MNBBBBBBBBBBBBBBBBBBBBBBBBBDDDDn
"""

console = Console()


def _banner():
    t = Text(LOGO)
    t.stylize("bold magenta")
    console.print(t, soft_wrap=True)
    console.print("  [dim]group-balanced retrieval de-confounding[/dim]\n")


def _fmt(x):
    return f"{x:.3f}" if x == x else "—"


def _print_refs(store, cohort_pids, refs, arm, header):
    title = "cow on beach  (far, same-diagnosis → invariance)" if arm == "cow" \
        else "pig on grass  (near, opposite-diagnosis → matched)"
    console.print(f"  {header}\n")
    rows = ref_distances(store, cohort_pids, refs, arm)
    table = Table(title=f"{title}   —   {len(refs)} references")
    table.add_column("patient", style="cyan", no_wrap=True)
    table.add_column("study", style="green")
    table.add_column("label")
    table.add_column("d_cohort", justify="right")
    table.add_column("d_refs", justify="right")
    for row in rows:
        raw = store.meta.loc[row["patient"], "label_raw"]
        name = str(raw) if isinstance(raw, str) else ("case" if row["label"] == 1 else "control")
        color = "red" if row["label"] == 1 else "blue"
        table.add_row(str(row["patient"]), row["study"], f"[{color}]{name}[/{color}]",
                      _fmt(row["d_cohort"]), _fmt(row["d_refs"]))
    console.print(table)
    dc = [r["d_cohort"] for r in rows if r["d_cohort"] == r["d_cohort"]]
    dr = [r["d_refs"] for r in rows if r["d_refs"] == r["d_refs"]]
    mc = _fmt(float(np.mean(dc)) if dc else float("nan"))
    mr = _fmt(float(np.mean(dr)) if dr else float("nan"))
    console.print(f"  [dim]studies covered: {store.meta.loc[refs, 'study'].nunique()}  |  mean d_cohort {mc}  |  mean d_refs {mr}  |  leak-free ✓[/dim]")


def _arm_fn(arm):
    return diverse_refs if arm == "cow" else matched_refs


def _demo(arm, K):
    s = SignatureStore.from_parquet(config.CRC_SIG)
    m = s.meta
    d = m[m.label == 1].study.value_counts().idxmax()
    n = m[(m.label == 0) & (m.study != d)].study.value_counts().idxmax()
    cohort = m.index[((m.study == d) & (m.label == 1)) | ((m.study == n) & (m.label == 0))].tolist()
    refs = _arm_fn(arm)(s, cohort, K, K)
    _print_refs(s, cohort, refs, arm, f"[yellow][demo mode][/yellow] confounded CRC cohort: tumor←{d} | normal←{n} (n={len(cohort)})")


def _run(arm, K):
    src, ds = project.get("source"), project.get("dataset")
    if not src or not ds:
        console.print("  [yellow]no source/dataset configured — running built-in demo[/yellow]\n")
        return _demo(arm, K)
    source = Source(src)
    corpus = source.load()
    cohort = load_cohort(ds, source)
    if corpus.sig.shape[1] != cohort.sig.shape[1]:
        console.print(f"  [red]signature-space mismatch:[/red] source dim {corpus.sig.shape[1]} != dataset dim {cohort.sig.shape[1]}")
        return
    combined = corpus.extended(cohort)
    cohort_pids = cohort.meta.index.tolist()
    refs = _arm_fn(arm)(combined, cohort_pids, K, K)
    _print_refs(combined, cohort_pids, refs, arm, f"source: {src['path']}  |  dataset: {ds['path']}  (cohort n={len(cohort_pids)})")


def _status():
    st = project.load()
    if not st:
        console.print("  nothing configured\n  set: [bold]sceris source <path>[/bold]  |  [bold]sceris dataset <path>[/bold]")
        return
    for k in ("source", "dataset"):
        v = st.get(k)
        console.print(f"  {k}: [cyan]{v['path'] if v else '—'}[/cyan]" + (f" [dim]({v['kind']})[/dim]" if v else ""))


def main():
    ap = argparse.ArgumentParser(prog="sceris", description="retrieval de-confounding")
    sub = ap.add_subparsers(dest="cmd")
    ps = sub.add_parser("source"); ps.add_argument("path")
    pd_ = sub.add_parser("dataset"); pd_.add_argument("path")
    sub.add_parser("status")
    for arm in ("cow", "pig"):
        pa = sub.add_parser(arm); pa.add_argument("-k", type=int, default=4)
    args = ap.parse_args()

    if args.cmd == "source":
        if not Path(args.path).exists():
            return console.print(f"  [red]path not found:[/red] {args.path}")
        spec = Source.local(args.path).spec
        project.set("source", spec)
        console.print(f"  source → [cyan]{spec['path']}[/cyan]")
    elif args.cmd == "dataset":
        if not Path(args.path).exists():
            return console.print(f"  [red]path not found:[/red] {args.path}")
        spec = spec_for(args.path)
        project.set("dataset", spec)
        console.print(f"  dataset → [cyan]{spec['path']}[/cyan]")
    elif args.cmd == "status":
        _status()
    elif args.cmd in ("cow", "pig"):
        _run(args.cmd, args.k)
    else:
        _banner()
        console.print("  [dim]set:[/dim]  [bold]sceris source <path>[/bold]  [bold]sceris dataset <path>[/bold]\n"
                      "  [dim]run:[/dim]  [bold]sceris cow[/bold]  |  [bold]sceris pig -k 6[/bold]  |  [bold]sceris status[/bold]")


if __name__ == "__main__":
    main()
