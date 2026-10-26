# Top-100 高置信未检出蛋白（假阳性候选）分析

- 数据：20 个共享种子的 `p_gated2`（MKMP-CLF）概率聚合。
- 筛选：y=0（20 种子均未检出）的蛋白中，按 mean(p_gated2) 降序取前 100。
- 校准质量：1131 个未检出蛋白中，仅 6 个 p>=0.9（0.5%）。

- 共 100 个蛋白，其中 99 个有 GO 注释。

## 校准质量（未检出蛋白的概率分布）

| 阈值 | 数量 | 占未检出比例 |
|---:|---:|---:|
| 0.99 | 2 | 0.2% |
| 0.95 | 3 | 0.3% |
| 0.9 | 6 | 0.5% |
| 0.85 | 9 | 0.8% |
| 0.8 | 15 | 1.3% |
| 0.7 | 27 | 2.4% |
| 0.5 | 50 | 4.4% |

## Top-100 假阳性候选列表

| accession | gene | p_gated2 | p_mkmp | p_esm | GO terms |
|---|---:|---:|---:|---:|---:|
| Q03393 | PTS | 0.9976 | 0.9934 | 0.9892 | 9 |
| P38606 | ATP6V1A | 0.9907 | 0.9592 | 0.9899 | 41 |
| P62847 | RPS24 | 0.9871 | 0.9313 | 0.9941 | 20 |
| Q9NSA1 | FGF21 | 0.9366 | 0.9499 | 0.5869 | 16 |
| P54289 | CACNA2D1 | 0.9104 | 0.9315 | 0.7499 | 26 |
| A6ND91 | ASPDH | 0.9017 | 0.8070 | 0.9239 | 4 |
| P10606 | COX5B | 0.8964 | 0.9617 | 0.3983 | 12 |
| P0DPB6 | POLR1D | 0.8543 | 0.7415 | 0.7990 | 17 |
| Q8IWX7 | UNC45B | 0.8517 | 0.8566 | 0.7548 | 9 |
| P17752 | TPH1 | 0.8334 | 0.7188 | 0.8526 | 12 |
| Q9BTY7 | HGH1 | 0.8328 | 0.6901 | 0.8941 | 2 |
| P13726 | F3 | 0.8247 | 0.7229 | 0.7186 | 27 |
| Q9BWT6 | MND1 | 0.8155 | 0.7134 | 0.8619 | 5 |
| Q13724 | MOGS | 0.8116 | 0.7082 | 0.8178 | 11 |
| A8MPS7 | YDJC | 0.8101 | 0.7736 | 0.7606 | 3 |
| P78352 | DLG4 | 0.7860 | 0.7958 | 0.4872 | 81 |
| O14977 | AZIN1 | 0.7810 | 0.8552 | 0.5259 | 10 |
| P18621 | RPL17 | 0.7646 | 0.6361 | 0.8361 | 16 |
| P31415 | CASQ1 | 0.7482 | 0.9213 | 0.3550 | 29 |
| Q5VU65 | NUP210L | 0.7417 | 0.7614 | 0.6399 | 2 |
| P22607 | FGFR3 | 0.7350 | 0.5434 | 0.6827 | 38 |
| Q13093 | PLA2G7 | 0.7328 | 0.7136 | 0.5351 | 18 |
| Q8TBF2 | PRXL2B | 0.7292 | 0.6608 | 0.7379 | 12 |
| O15078 | CEP290 | 0.7206 | 0.6444 | 0.5171 | 31 |
| Q96SB8 | SMC6 | 0.7140 | 0.7491 | 0.5762 | 27 |
| Q08AG7 | MZT1 | 0.7133 | 0.5770 | 0.8435 | 11 |
| O75794 | CDC123 | 0.7013 | 0.6764 | 0.5714 | 7 |
| O60499 | STX10 | 0.6903 | 0.6036 | 0.7037 | 20 |
| Q9H3J6 | MTRFR | 0.6846 | 0.6690 | 0.5801 | 10 |
| P53384 | NUBP1 | 0.6613 | 0.7418 | 0.5716 | 18 |
| A6NFN3 | RBFOX3 | 0.6555 | 0.8072 | 0.2574 | 9 |
| P07101 | TH | 0.6526 | 0.6739 | 0.4135 | 42 |
| Q5VV43 | KIAA0319 | 0.6524 | 0.4335 | 0.6373 | 11 |
| P23458 | JAK1 | 0.6403 | 0.5911 | 0.6168 | 52 |
| P16671 | CD36 | 0.6374 | 0.7979 | 0.4326 | 106 |
| Q14690 | PDCD11 | 0.6348 | 0.4892 | 0.6308 | 13 |
| O94906 | PRPF6 | 0.6177 | 0.7346 | 0.3866 | 27 |
| Q9Y5V0 | ZNF706 | 0.5949 | 0.7921 | 0.3660 | 5 |
| O00294 | TULP1 | 0.5936 | 0.5666 | 0.4047 | 19 |
| O60869 | EDF1 | 0.5857 | 0.4905 | 0.6811 | 15 |
| Q8TDH9 | BLOC1S5 | 0.5762 | 0.5647 | 0.4905 | 17 |
| O60870 | KIN | 0.5699 | 0.4803 | 0.5777 | 17 |
| Q12794 | HYAL1 | 0.5613 | 0.6309 | 0.5148 | 28 |
| P52815 | MRPL12 | 0.5564 | 0.3992 | 0.5804 | 13 |
| Q8N960 | CEP120 | 0.5349 | 0.3801 | 0.5716 | 12 |
| O76093 | FGF18 | 0.5323 | 0.4894 | 0.5565 | 18 |
| Q9HAV7 | GRPEL1 | 0.5265 | 0.5163 | 0.4276 | 13 |
| Q9UJA5 | TRMT6 | 0.5258 | 0.4946 | 0.3749 | 7 |
| P04844 | RPN2 | 0.5078 | 0.3018 | 0.7969 | 10 |
| Q01968 | OCRL | 0.5026 | 0.5612 | 0.3689 | 41 |
| P46736 | BRCC3 | 0.5000 | 0.4294 | 0.4646 | 32 |
| O75330 | HMMR | 0.4995 | 0.3877 | 0.5422 | 14 |
| Q9Y285 | FARSA | 0.4979 | 0.3562 | 0.7439 | 15 |
| P20648 | ATP4A | 0.4922 | 0.5023 | 0.5335 | 18 |
| Q13601 | KRR1 | 0.4919 | 0.2506 | 0.7757 | 12 |
| Q9NX04 | AIRIM | 0.4806 | 0.6361 | 0.3372 | 10 |
| Q14206 | RCAN2 | 0.4764 | 0.4703 | 0.4841 | 5 |
| Q709F0 | ACAD11 | 0.4745 | 0.5552 | 0.4656 | 17 |
| Q7Z3B4 | NUP54 | 0.4675 | 0.5316 | 0.3149 | 10 |
| Q08AF3 | SLFN5 | 0.4623 | 0.7222 | 0.2438 | 4 |
| P29475 | NOS1 | 0.4566 | 0.3537 | 0.4937 | 78 |
| Q9Y2G2 | CARD8 | 0.4565 | 0.4149 | 0.3224 | 31 |
| Q6RSH7 | VHLL | 0.4533 | 0.5942 | 0.2077 | 1 |
| Q6UWF7 | NXPE4 | 0.4440 | 0.2355 | 0.7339 | 2 |
| P27986 | PIK3R1 | 0.4356 | 0.5533 | 0.3458 | 63 |
| P57060 | RWDD2B | 0.4339 | 0.5998 | 0.2514 | 1 |
| O95470 | SGPL1 | 0.4311 | 0.5045 | 0.4369 | 11 |
| Q96C36 | PYCR2 | 0.4200 | 0.4618 | 0.3777 | 10 |
| P03971 | AMH | 0.4174 | 0.5418 | 0.2616 | 22 |
| P35916 | FLT4 | 0.4161 | 0.3333 | 0.4438 | 38 |
| P21781 | FGF7 | 0.4131 | 0.3033 | 0.4571 | 30 |
| Q6BDS2 | BLTP3A | 0.4111 | 0.4117 | 0.3030 | 4 |
| Q96K12 | FAR2 | 0.4105 | 0.5974 | 0.2579 | 7 |
| P19235 | EPOR | 0.4096 | 0.5252 | 0.2128 | 12 |
| Q9BQW3 | EBF4 | 0.4079 | 0.3538 | 0.3430 | 14 |
| O75694 | NUP155 | 0.4073 | 0.3663 | 0.4823 | 14 |
| P20783 | NTF3 | 0.4070 | 0.3247 | 0.4484 | 24 |
| P51116 | FXR2 | 0.4067 | 0.1037 | 0.6397 | 28 |
| Q9BZI7 | UPF3B | 0.4063 | 0.3290 | 0.4246 | 16 |
| Q7Z392 | TRAPPC11 | 0.4062 | 0.4177 | 0.5172 | 11 |
| Q6BDI9 | REP15 | 0.4057 | 0.4181 | 0.3824 | 7 |
| Q16595 | FXN | 0.4052 | 0.5087 | 0.2709 | 27 |
| Q9BW85 | YJU2 | 0.4043 | 0.3113 | 0.4653 | 9 |
| Q9BTT4 | MED10 | 0.3959 | 0.4546 | 0.4607 | 13 |
| A8MT69 | CENPX | 0.3953 | 0.2187 | 0.6974 | 18 |
| Q8NBM8 | PCYOX1L | 0.3780 | 0.2981 | 0.5894 | 9 |
| Q96CV9 | OPTN | 0.3772 | 0.3809 | 0.3275 | 30 |
| P48357 | LEPR | 0.3728 | 0.3158 | 0.3588 | 37 |
| Q08379 | GOLGA2 | 0.3726 | 0.3871 | 0.3201 | 31 |
| Q14498 | RBM39 | 0.3670 | 0.1702 | 0.5538 | 14 |
| Q8N8E3 | CEP112 | 0.3666 | 0.3205 | 0.4093 | 10 |
| Q07020 | RPL18 | 0.3625 | 0.1817 | 0.5905 | 18 |
| Q13114 | TRAF3 | 0.3602 | 0.4430 | 0.3652 | 40 |
| P38484 | IFNGR2 | 0.3575 | 0.3695 | 0.3974 | 16 |
| P32971 | TNFSF8 | 0.3521 | 0.4812 | 0.2495 | 11 |
| Q9BQB4 | SOST | 0.3515 | 0.3520 | 0.4210 | 17 |
| P35913 | PDE6B | 0.3365 | 0.6039 | 0.1922 | 19 |
| O14920 | IKBKB | 0.3292 | 0.3853 | 0.3596 | 47 |
| P0CW23 | AKAIN1 | 0.3282 | 0.4034 | 0.3497 | 4 |
| Q9ULN7 | PNMA8B | 0.3252 | 0.5190 | 0.1705 | 0 |

## GO 富集（FDR<0.05）

| GO | term | namespace | fg | bg | FDR |
|---|---|---|---:|---:|---:|
| GO:0042462 | eye photoreceptor cell development | biological_process | 3 | 10 | 3.31e-03 |
| GO:0036228 | protein localization to nuclear inner membrane | biological_process | 2 | 4 | 5.70e-03 |
| GO:0004896 | cytokine receptor activity | molecular_function | 4 | 50 | 6.66e-03 |
| GO:0016714 | oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced pteridine as one donor, and incorporation of one atom of oxygen | molecular_function | 2 | 4 | 6.84e-03 |
| GO:0008543 | fibroblast growth factor receptor signaling pathway | biological_process | 4 | 43 | 7.33e-03 |
| GO:0019903 | protein phosphatase binding | molecular_function | 5 | 92 | 7.73e-03 |
| GO:0043410 | positive regulation of MAPK cascade | biological_process | 6 | 178 | 8.06e-03 |
| GO:0005111 | type 2 fibroblast growth factor receptor binding | molecular_function | 2 | 5 | 8.12e-03 |
| GO:0008083 | growth factor activity | molecular_function | 5 | 139 | 1.43e-02 |
| GO:0005515 | protein binding | molecular_function | 82 | 13384 | 1.52e-02 |
| GO:0019221 | cytokine-mediated signaling pathway | biological_process | 5 | 136 | 1.58e-02 |
| GO:1902514 | regulation of calcium ion transmembrane transport via high voltage-gated calcium channel | biological_process | 2 | 10 | 1.80e-02 |
| GO:0090161 | Golgi ribbon formation | biological_process | 2 | 11 | 1.80e-02 |
| GO:0035666 | TRIF-dependent toll-like receptor signaling pathway | biological_process | 2 | 11 | 1.92e-02 |
| GO:0034121 | regulation of toll-like receptor signaling pathway | biological_process | 2 | 10 | 1.93e-02 |
| GO:0035631 | CD40 receptor complex | cellular_component | 2 | 11 | 2.04e-02 |
| GO:0009072 | aromatic amino acid metabolic process | biological_process | 2 | 10 | 2.10e-02 |
| GO:0030315 | T-tubule | cellular_component | 3 | 49 | 2.20e-02 |
| GO:0060333 | type II interferon-mediated signaling pathway | biological_process | 2 | 13 | 2.27e-02 |
| GO:0038060 | nitric oxide-cGMP-mediated signaling | biological_process | 2 | 13 | 2.40e-02 |
| GO:0016529 | sarcoplasmic reticulum | cellular_component | 3 | 53 | 2.63e-02 |
| GO:0031410 | cytoplasmic vesicle | cellular_component | 7 | 369 | 2.73e-02 |
| GO:1990229 | iron-sulfur cluster assembly complex | cellular_component | 2 | 16 | 2.74e-02 |
| GO:0046677 | response to antibiotic | biological_process | 2 | 16 | 2.86e-02 |
| GO:0030214 | hyaluronan catabolic process | biological_process | 2 | 17 | 2.97e-02 |
| GO:0003723 | RNA binding | molecular_function | 15 | 1401 | 3.47e-02 |
| GO:0006941 | striated muscle contraction | biological_process | 3 | 63 | 3.47e-02 |
| GO:0005790 | smooth endoplasmic reticulum | cellular_component | 2 | 20 | 3.54e-02 |
| GO:0005829 | cytosol | cellular_component | 38 | 5095 | 3.55e-02 |
| GO:0007267 | cell-cell signaling | biological_process | 5 | 221 | 3.97e-02 |
