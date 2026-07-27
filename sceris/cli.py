import argparse
from rich.console import Console
from rich.table import Table
from rich.text import Text
from . import config
from .store import SignatureStore
from .retrieval import diverse_refs, matched_refs

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


def _demo(arm, K):
    console.print(f"[dim]loading CRC signatures…[/dim]")
    s = SignatureStore.from_parquet(config.CRC_SIG)
    m = s.meta
    d = m[m.label == 1].study.value_counts().idxmax()
    n = m[(m.label == 0) & (m.study != d)].study.value_counts().idxmax()
    cohort = m.index[((m.study == d) & (m.label == 1)) | ((m.study == n) & (m.label == 0))].tolist()
    console.print(f"  confounded cohort: [red]tumor[/red]←{d}  |  [blue]normal[/blue]←{n}  (n={len(cohort)})\n")

    fn = diverse_refs if arm == "cow" else matched_refs
    refs = fn(s, cohort, K, K)
    title = "🐄  cow on beach  (far, same-diagnosis → invariance)" if arm == "cow" \
        else "🐷  pig on grass  (near, opposite-diagnosis → matched)"
    table = Table(title=f"{title}   —   {len(refs)} references")
    table.add_column("patient", style="cyan", no_wrap=True)
    table.add_column("study", style="green")
    table.add_column("label")
    for r in refs:
        row = s.meta.loc[r]
        lab = "[red]tumor[/red]" if row.label == 1 else "[blue]normal[/blue]"
        table.add_row(str(r), str(row.study), lab)
    console.print(table)
    console.print(f"  [dim]studies covered: {s.meta.loc[refs, 'study'].nunique()}  |  leak-free ✓[/dim]")


def main():
    ap = argparse.ArgumentParser(prog="sceris", description="retrieval de-confounding, in the terminal")
    ap.add_argument("arm", nargs="?", choices=["cow", "pig"], help="which retrieval arm to demo")
    ap.add_argument("-k", type=int, default=4, help="references per class (default 4)")
    args = ap.parse_args()
    _banner()
    if args.arm:
        _demo(args.arm, args.k)
    else:
        console.print("  [dim]try:[/dim]  [bold]sceris cow[/bold]   |   [bold]sceris pig -k 6[/bold]")


if __name__ == "__main__":
    main()
