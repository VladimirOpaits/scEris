from rich.console import Console
from rich.text import Text

console = Console()

logo = r"""
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

text = Text(logo)
text.stylize("bold magenta")
console.print(text)