### Hybrid bond number

| case number  | Triton | Simple |
| ------------ | ------ | ------ |
| case2        | 149    | 332    |
| case2 hidden | 149    | 332    |
| case3        | 1070   | 9187   |
| case3 hidden | 1086   | 9187   |
| case4        | 26602  | 5589   |
| case4 hidden | 26764  | 5589   |

Set pins for Input Output info by fanout number



### Score

Env: default
density: 1.0, no net weight

| case number  | Beta 1st  | Final 1st | Mine - initial | Mine Final |
| ------------ | --------- | --------- | -------------- | ---------- |
| case2        | 2189151   | 1960913   | 1901792        | 2214728    |
| case2 hidden | 2958953   | 2555461   | 2220967        | 2596823    |
| case3        | 36191154  | 30247740  | 30837439       | 35413604   |
| case3 hidden | 31911473  | 27650329  | 27734508       | 31396814   |
| case4        | 328483651 | 274026678 | 317682092      | 363167567  |
| case4 hidden | 361787452 | 301193374 | 331920892      | 382798200  |

Change Env:
Target density from 1.0 to 1.5, net weight 1.5 for intersected net
This follows the commit which has the parent commit of `14380f90`

| case number  | Final 1st | Mine Final |
| ------------ | --------- | ---------- |
| case2        | 1960913   | 1846984    |
| case2 hidden | 2555461   | 2172458    |
| case3        | 30247740  | 28780432   |
| case3 hidden | 27650329  | 25646473   |
| case4        | 274026678 | 301201197  |
| case4 hidden | 301193374 | 323151813  |









### SCREEN NUMBER

| screen ID                  | Case number   |
| -------------------------- | ------------- |
| 2855955.pts-54.csdl-ubuntu | case 2        |
| 4125299.pts-54.csdl-ubuntu | case 2 hidden |
| 21616.pts-54.csdl-ubuntu   | case 3        |
| 555901.pts-54.csdl-ubuntu  | case 3 hidden |
| 1801352.pts-54.csdl-ubuntu | case 4        |
| 1871183.pts-54.csdl-ubuntu | case 4 hidden |