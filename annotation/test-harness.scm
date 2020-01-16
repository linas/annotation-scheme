;
: Quick and dirty test harness for exploring computational
; performance of the various different processing pathways.
;
; This is the code that was used to obtain the results reported
; in https://github.com/MOZI-AI/annotation-scheme/issues/98
;


; Basic list of 681 genes
; See also lmpd-genes for a list of 1482
(define gene-list (list "TSPAN6" "NDUFAF7" "RBM5" "SLC7A2" "NDUFAB1"
"DVL2" "SKAP2" "DHX33" "MSL3" "BZRAP1" "GTF2IRD1" "IL32" "RPS20" "SCMH1"
"CLCN6" "RNF14" "ATP2C1" "IGF1" "GLRX2" "FAS" "ATP6V0A1" "FBXO42" "JADE2"
"PREX2" "NOP16" "LMO3" "R3HDM1" "ERCC8" "HOMER3" "USE1" "OPN3" "SZRD1"
"ATG5" "CAMK2B" "MPC1" "MRPS24" "ZNF275" "TAF2" "TAF11" "IPO5" "NDUFB4"
"DIP2B" "MPPED2" "IARS2" "ERLEC1" "UFD1L" "PDCD2" "ACADVL" "ENO1" "FRYL"
"SEC31B" "KIFAP3" "NT5C2" "GPC4" "ITGA8" "PPP2R5C" "RBFOX1" "ITM2A"
"NRD1" "VDAC3" "CBFA2T2" "FKBP7" "SAR1A" "DUSP13" "PGR" "EPB41L3" "OXCT1"
"SLC27A5" "WBP11" "NCOA1" "MAPRE3" "MGST2" "DIMT1" "RBM22" "TMED2"
"HUWE1" "NLK" "UIMC1" "GNAS" "COQ9" "NSFL1C" "TASP1" "MRPS33" "NDUFB2"
"TXNL1" "MYL6" "HDAC6" "DHPS" "CREM" "PSMD8" "CIRBP" "HNRNPM" "SF3A1"
"POLR2F" "HMGXB4" "CHKB" "ZMAT5" "RBM23" "VTI1B" "TIMM9" "GSTZ1"
"RPS6KA5" "PSMB5" "PFDN4" "PSMA7" "NDUFAF5" "ATP1B4" "ALG13" "SUV39H1"
"SCML2" "PGK1" "KLF5" "TSC22D1" "MGRN1" "SLC7A6" "CMC2" "CPPED1" "MYEF2"
"CPQ" "LEPROTL1" "PPP2CB" "KLHDC4" "POP4" "AKT2" "RABAC1" "CARD8" "PON2"
"SSBP1" "BUD31" "MEST" "CHCHD3" "COA1" "BLVRA" "PLGRKT" "BAG1" "EXOSC3"
"RASSF4" "KAZALD1" "PITRM1" "EBF3" "LGI1" "MTMR4" "CDK5RAP3" "ENO3"
"ICAM2" "EZH1" "MRPL27" "HDAC5" "DUSP3" "DCUN1D4" "NDUFC1" "ZBTB16"
"COMMD9" "ATP5B" "ELK3" "ALDH2" "STX2" "GPR133" "MRPL51" "GAPDH" "TPI1"
"TMEM14C" "GCNT2" "NCOA7" "FANCE" "E2F3" "ACOT13" "COX7A2" "ENPP5"
"PCDHB2" "TMCO6" "TTC1" "POLR3G" "BNIP1" "TIMMDC1" "BCL6" "FAM162A"
"PRKAR2A" "TUSC2" "SSR3" "MOB1A" "NCL" "MPV17" "HSPE1" "ZNF142" "ID2"
"PNO1" "TMEM59" "LAMTOR2" "LEPR" "CTH" "TMEM9" "MRPS15" "SDHB" "FAAH"
"DPH5" "ANKRD13C" "VAMP8" "NDUFB3" "GDA" "MAPKAP1" "YPEL5" "RCL1" "GRIA2"
"PCMT1" "KIAA1217" "PAIP2" "ARL1" "SOCS2" "ECHDC2" "A1BG" "ZNF211"
"GTDC1" "CCDC18" "HNRNPA2B1" "FAM126A" "CLTA" "CISD1" "CDKN2C" "RASSF8"
"ATP7B" "ITIH5" "SMUG1" "NDUFAF4" "PLA2G12A" "PFKFB2" "ATP5E" "STAMBP"
"SNRNP27" "NQO2" "EMC3" "TMTC4" "SNRPD2" "ID1" "KDM5C" "ST3GAL3" "EMC1"
"UQCR11" "RNF6" "SHFM1" "STEAP4" "RBM48" "TUBGCP6" "EMC4" "RABL5" "MTX2"
"TXNDC17" "DAD1" "ECSIT" "CDC16" "RPL36" "MRPL34" "LSM4" "CYP2E1"
"ZNF337" "PRRC2B" "COX4I1" "NFATC1" "PDLIM4" "PSME3" "NDUFA2" "RAF1"
"ENOSF1" "MRPL35" "SERPINF1" "GRSF1" "WBP2" "PRMT7" "POMP" "MYH8" "BEX2"
"ACTR3B" "ARL8B" "EMC7" "LAMTOR5" "KCTD1" "DDB2" "PUM1" "NREP" "CTSL"
"FAM189A2" "MDFIC" "CAPRIN1" "CD63" "TSPAN31" "MRPL44" "NDUFB5" "TANK"
"VPS45" "PSMB7" "UBAP2" "GRHPR" "BPHL" "UQCC2" "NUMA1" "DCUN1D5" "UNC13C"
"TUBGCP4" "SLC28A2" "CYP1B1" "COX17" "PARP16" "HERC5" "PPA2" "LGR5"
"NEDD1" "SLAIN1" "GRTP1" "TMX1" "SERF2" "SRP14" "WDR61" "TPM1" "SEC11A"
"PMM2" "MYLK3" "CLTC" "GAREM" "GREB1L" "RNF165" "TXNL4A" "NFIC" "PSMB6"
"PSMA5" "CELSR2" "TIPRL" "UFC1" "SETDB1" "ADAMTSL4" "JTB" "HAX1" "ACP1"
"NVL" "DEGS1" "PSEN2" "PDIA6" "SFXN5" "ZNF385B" "UBR3" "GULP1" "ATG3"
"SRPRB" "TMEM108" "SLIT2" "DDIT4L" "PAM" "TSLP" "ATG12" "BOD1" "MYLK4"
"DYNLT1" "PSPH" "ATXN7L1" "ZMYM3" "ARHGAP36" "HMBOX1" "PXDNL" "GOLGA7"
"POLR2K" "EIF3H" "NDUFB9" "TATDN1" "FAM171A1" "FAM188A" "GSTO1" "EIF3M"
"ARFGAP2" "DAK" "CPSF7" "CABLES2" "SAP18" "HNMT" "TIMM8B" "PTS" "NDUFC2"
"BICD1" "ZNF385D" "GEMIN6" "ATP5A1" "CCDC50" "UTRN" "ZNF117" "RASSF3"
"HNRNPU" "TRIP12" "LGI4" "MSI2" "UCHL1" "ATP5G3" "TCEB1" "PSMA8" "DPH3"
"ACSS1" "TMEM55A" "GOLGA7B" "CNOT8" "NUP205" "KCNMA1" "SCAF4" "MALSU1"
"PDE6D" "MUM1L1" "AFAP1L1" "WDR19" "MRPL17" "CNNM4" "ZFAND2B" "AGPAT6"
"SNF8" "CBR1" "TPPP3" "CCDC28B" "SSU72" "PCNT" "PCYT1A" "COX7A1" "BCL6B"
"POLR3K" "NXF1" "OLFML2B" "MRPL55" "RFTN2" "VSNL1" "RPRD2" "BOLA3"
"MRPS18C" "ARPC2" "SUCLG1" "PPM1L" "ADAMTS9" "ELP6" "SMIM12" "SFMBT1"
"RAD54L2" "HPGD" "ARFIP1" "UQCRQ" "HCN1" "IQUB" "SUN1" "DCAF13" "NDUFB6"
"CFL2" "KLHDC2" "TRUB1" "ZFYVE1" "PSMC3" "DPCD" "ALKBH3" "CCT2" "ZNF202"
"USP54" "NDST2" "SEC11C" "CENPV" "HSP90B1" "AP1G1" "PPIB" "FAM96A"
"SMAD3" "PDIA3" "FBXO22" "ATP5L" "GPX4" "ATP5H" "POLR2G" "COPS6" "RAB4A"
"DNAJC7" "COA6" "MMADHC" "MPLKIP" "HEXIM2" "COL3A1" "CNNM3" "USP39"
"RNF181" "LRRC28" "PPIC" "STXBP6" "MFF" "ATP5I" "PARM1" "NSMCE1" "EFNA1"
"CRADD" "SIN3A" "CNBP" "ZNF32" "TOR1AIP2" "BRD3" "SF3B5" "STX8" "UBB"
"DNAJC18" "NUDCD2" "GPR27" "KBTBD2" "TRIAP1" "MTSS1" "ZNF160" "FEZ2"
"FAM86JP" "RBKS" "EXOSC1" "PDE7B" "MRPL36" "BPTF" "ZNF540" "EXOSC10"
"FRMD3" "TP53RK" "SYNPO2" "MYEOV2" "MRPL52" "TMEM134" "BNC2" "KDM2A"
"RHOD" "COMMD1" "GLRX" "DMRT2" "FAM86B3P" "MINOS1" "UQCRH" "PHC3" "USMG5"
"HOXB2" "NRROS" "MSRB3" "SNAPC5" "MRPL11" "IL20RB" "AKIRIN1" "TMEM167A"
"NDUFA11" "CADM2" "LSM1" "TMEM9B" "SCUBE2" "UCP3" "PAAF1" "MRPL48" "RTTN"
"KCMF1" "RMDN1" "IRX5" "MAMSTR" "POLE" "FAM210A" "ATOX1" "KIAA0195"
"MLF1" "GRAMD1C" "BOLA1" "TMEM11" "RIIAD1" "SEPW1" "MAGED1" "FCER1A"
"PSMG4" "SSR4" "TRAPPC5" "PLAG1" "LSM10" "COA4" "NEUROG1" "MRPS11"
"MRPS16" "FAM104B" "SATB1" "NDN" "CNOT10" "ZNF662" "CEP57L1" "SFXN4"
"FAM162B" "DENND5A" "UBE2F" "PCDH9" "MROH7" "NDUFA12" "APOO" "PTRHD1"
"RPS27L" "ADSSL1" "UBALD2" "ROR1" "NR2F2" "PSMD13" "ANKFY1" "ZNF529"
"POLR1D" "TOR3A" "PPP1CC" "RGS9BP" "KRT10" "GPATCH8" "RPS19BP1" "CMC1"
"MAGI2" "EXD3" "MORN2" "COL4A5" "COMMD6" "S100A16" "RNFT1" "HN1"
"BLOC1S2" "LCOR" "XPNPEP3" "ACN9" "TECPR2" "TOMM7" "ADA" "NOP9" "CASP4"
"METTL9" "ZNF398" "ZNF682" "MRPL21" "HTT" "COL4A6" "PHF2" "FAM118B"
"FAM49A" "MRPL42" "UBL5" "CCDC69" "NCOA6" "MT-ND5" "ZNF521" "MT-ND3"
"SELT" "OSTC" "TSEN15" "MT-ND4" "DCLRE1A" "CCDC167" "DMD" "SNORA63"
"TATDN3" "FAM229B" "INPP5B" "COL15A1" "ZBTB48" "ZNF783" "FLJ00104"
"SARNP" "DNAJC19" "SNORA12" "U3" "U3" "SUPT4H1" "ANKRD39" "NDUFS3"
"SLC23A3" "ARL16" "ZSWIM7" "COL28A1" "SLC35E2" "TENM3" "FAM19A5" "TPI1P1"
"OST4" "EEF1DP3" "HSBP1" "MCTS1" "DICER1-AS1" "CDKN2AIPNL" "FAM200B"
"GOLGA2P5" "MRPL33" "NME1-NME2" "UBE2V1" "SMIM20" "APOPT1" "CUX1"
"UBE2F-SCLY" "GOLGA7B" "NCBP2-AS2"))

(use-modules (annotation biogrid))
(use-modules (annotation gene-go))
(use-modules (annotation gene-pathway))
(use-modules (annotation functions))
(use-modules (annotation util))

; (load "inst.scm")
; (define smpdb-ctr (accum-time "smpdb"))
; (define reactome-ctr (accum-time "reactome"))
; (define find-pathway-genes-ctr (accum-time "find-pathway-genes"))
; (define add-pathway-genes-ctr (accum-time "add-pathway-genes"))
; (define find-go-term-ctr (accum-time "find-go-term"))
; (define find-memberln-ctr (accum-time "find-memberln"))
; (define add-go-info-ctr (accum-time "add-go-info"))
; (define find-parent-ctr (accum-time "find-parent"))

(define (report)
;pathway stuff
(smpdb-ctr #:report? #t)
(reactome-ctr #:report? #t)
(find-pathway-genes-ctr #:report? #t)
(add-pathway-genes-ctr #:report? #t)
(find-pathway-member-ctr #:report? #t)
(pathway-gene-interactors-ctr #:report? #t)
(generate-interactors-ctr #:report? #t)
(pathway-hierarchy-ctr #:report? #t)
(check-pathway-ctr #:report? #t)
(find-protein-ctr #:report? #t)
(find-mol-ctr #:report? #t)
(find-go-term-ctr #:report? #t)

; grid stuff
(match-gene-interactors-ctr #:report? #t)
(find-output-interactors-ctr #:report? #t)
(generate-result-ctr #:report? #t)
(build-interaction-ctr #:report? #t)
(find-protein-form-ctr #:report? #t)

; common to grid and path
(find-name-ctr #:report? #t)
(find-pubmed-id-ctr #:report? #t)

; ???
(find-memberln-ctr #:report? #t)
(add-go-info-ctr #:report? #t)
(find-parent-ctr #:report? #t)
(locate-node-ctr #:report? #t)
(add-loc-ctr #:report? #t)
)


(define (do-anno nparents)
	(define start (current-time))
	(define anno (gene-go-annotation gene-list
		"my-go-anno-results"
		#:parents nparents))
	(define elapse (- (current-time) start))
	(format #t "GO Annotation took ~A seconds\n" elapse)
	*unspecified*
)

(define (do-path nparents)
	(define start (current-time))
	(define anno (gene-pathway-annotation gene-list
		"path-results"
		#:pathway "reactome smpdb"
		; XXX if the pathway is set, then it crashes, see issue #91
		; #:namespace "biological_process molecular_function cellular_component"
		#:parents nparents))
	(define elapse (- (current-time) start))
	(format #t "Path Annotation (p=~A) took ~A seconds\n" nparents elapse)
	*unspecified*
)

(define (do-path-ns nparents)
	(define start (current-time))
	(define anno (gene-pathway-annotation gene-list
		"path-with-ns-results"
		#:pathway "reactome smpdb"
		#:namespace "biological_process molecular_function cellular_component"
		#:parents nparents))
	(define elapse (- (current-time) start))
	(format #t "Path Annotation (p=~A) took ~A seconds\n" nparents elapse)
	*unspecified*
)

(define (do-one-path gename)
	(define start (current-time))
	(define anno (gene-pathway-annotation (list gename)
		"do-one-path-results"
		#:pathway "reactome smpdb"
		; XXX if the pathway is set, then it crashes, see issue #91
		; #:namespace "biological_process molecular_function cellular_component"
		#:parents 0))
	(define elapse (- (current-time) start))
	(format #t "Path Annotation for ~A took ~A seconds; got ~A annotations\n"
		gename elapse (cog-arity anno))
	*unspecified*
)

(define (do-grid-protein nparents)
	(define start (current-time))
	(define anno (biogrid-interaction-annotation gene-list
		"my-biogrid-results"
		#:namespace "biological_process molecular_function cellular_component"
		#:interaction "Proteins"
		; #:interaction "Genes"
		#:parents nparents))
	(define elapse (- (current-time) start))
	(format #t "Grid Annotation took ~A seconds\n" elapse)
	*unspecified*
)

(define (do-grid-gene nparents)
	(define start (current-time))
	(define anno (biogrid-interaction-annotation gene-list
		"my-biogrid-results"
		#:namespace "biological_process molecular_function cellular_component"
		#:interaction "Genes"
		#:parents nparents))
	(define elapse (- (current-time) start))
	(format #t "Grid Annotation took ~A seconds\n" elapse)
	*unspecified*
)

; (start-cogserver "path.conf")
; (start-cogserver "grid.conf")