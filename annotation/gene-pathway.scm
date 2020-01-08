;;; MOZI-AI Annotation Scheme
;;; Copyright © 2019 Abdulrahman Semrie
;;; Copyright © 2019 Hedra Seid
;;;
;;; This file is part of MOZI-AI Annotation Scheme
;;;
;;; MOZI-AI Annotation Scheme is free software; you can redistribute
;;; it and/or modify it under the terms of the GNU General Public
;;; License as published by the Free Software Foundation; either
;;; version 3 of the License, or (at your option) any later version.
;;;
;;; This software is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;; General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this software.  If not, see
;;; <http://www.gnu.org/licenses/>.

(define-module (annotation gene-pathway)
      #:use-module (annotation functions)
      #:use-module (annotation util)
      #:use-module (opencog)
      #:use-module (opencog exec)
      #:use-module (opencog bioscience)
      #:use-module (annotation parser)
      #:export (gene-pathway-annotation)
)

(use-modules (ice-9 format))
(define (accum-time name)
   (let ((fname name)
         (elapsed 0)
         (calls 0)
         (start 0))
      (lambda* (#:key (enter? #f) (report? #f))
         (if report?
            (if (< 0 calls)
               (format #t 
                  "Time: ~9f secs. calls: ~A avg: ~8,1f usec/call for ~A\n"
                  (* 1.0e-9 elapsed)
                  calls
                  (/ (* 1.0e-3 elapsed) calls)
                  fname)
               (format #t "Zero calls to ~A\n" fname))
            (if enter?
               (set! start (get-internal-real-time))
               (begin
                  (set! elapsed (+ elapsed
                     (- (get-internal-real-time) start)))
                  (set! calls (+ calls 1)))))))
)


(define-public reactome-ctr (accum-time "reactome"))
(define-public smpdb-ctr (accum-time "smpdb"))


(define report-gc
	(let* ((stime 0)
			(sngc  0)
			(srun (get-internal-run-time))
			(swall (get-internal-real-time))
			(igtime (cdr (assoc 'gc-time-taken (gc-stats))))
			(irun srun)
			(iwall swall)
		)
	(lambda ()
		(define agc (gc-stats))
		(define time (cdr (assoc 'gc-time-taken agc)))
		(define run (get-internal-run-time))
		(define wall (get-internal-real-time))
		(define heap (cdr (assoc 'heap-size agc)))
		(define ngc  (cdr (assoc 'gc-times agc)))
		(define nugc  (- ngc sngc))
		(define gsecs (/ (- time stime) 1000000000.))
		(define gacc (/ (- time igtime) 1000000000.))
		(define psecs (/ (- run irun) 1000000000.))
		(define wsecs (/ (- wall iwall) 1000000000.))
		(define hrs (/ gsecs 3600.))
		(define mb (/ heap 1000000.))
		; cum proc time has gc time subtracted.
		(format #t "Heap: ~6f MB  #gc: ~A ~5f secs cum: proc: ~7f wall: ~7f\n"
			mb nugc gsecs (- psecs gacc) wsecs)
		(set! stime time)
		(set! sngc ngc)
		(set! srun run)
		(set! swall wall)
	)
))


(define* (gene-pathway-annotation gene_nodes file-name #:key (pathway "reactome") (include_prot "True") (include_sm "True") (namespace "") (parents 0)  (biogrid 1))
(format #t "GC start: ~A\n" (gc-stats))
    (let (
[gctr 0]
[numg (length gene_nodes)]
[result '()]
          [pwlst '()]
          [go (if (string=? namespace "") (ListLink)
                (ListLink (ConceptNode namespace) (Number parents)))])

    (for-each (lambda (gene)
(set! gctr (+ 1 gctr))
      (set! result (append result (node-info (GeneNode gene))))
      (for-each (lambda (pathw)
          (if (equal? pathw "smpdb")
(let ((start (get-internal-real-time)))
              (set! result (append result (smpdb gene include_prot include_sm go biogrid)))
(format #t "Did smpdb ~A of ~A for ~A result-len=~A time=~6f\n"
gctr numg gene (length result) (* 1.0e-9 (- (get-internal-real-time) start)))
(report-gc)
)
              )
          (if (equal? pathw "reactome")
              (begin
              (let* ((start (get-internal-real-time))
[res (reactome gene include_prot include_sm pwlst go biogrid)])
(format #t "Did reactome, now append ~A ~A\n" (length (car res)) (length (cdr
res)))
                (set! result (append result (car res)))
                (set! pwlst (append pwlst (cdr res)))
(format #t "Did reactome ~A of ~A for ~A result-len=~A pwlen=~A time=~6f\n"
gctr numg gene (length result) (length pwlst) (* 1.0e-9 (- (get-internal-real-time) start)))
(report-gc)
              )))
          )(string-split pathway #\ ))
    ) gene_nodes)

(format #t "done them all\n")
    (let (
      [res (ListLink (ConceptNode "gene-pathway-annotation") (ListLink result))]
    )
      (write-to-file res file-name "gene-pathway")
      res
    )
))


;; From SMPDB

(define (smpdb gene prot sm go biogrid)
(smpdb-ctr #:enter? #t)
  (let (
    [pw (find-pathway-member (GeneNode gene) "SMP")]
    [ls '()]
  )

  (set! ls (flatten (map (lambda (path)
      (let (
        [node (cog-outgoing-atom (cog-outgoing-atom path 0) 1)]
        [tmp '()]
      )
      (if (equal? sm "True")
          (set! tmp (append tmp (cog-outgoing-set (find-mol node "ChEBI"))))
      )
      (set! tmp (append tmp (find-pathway-genes node go)))
      (if (equal? prot "True")
        (let ([prots (cog-outgoing-set (find-mol node "Uniprot"))])
          (if (not (null? prots))
            (set! tmp (append tmp prots))
            (set! tmp (append tmp (node-info node))))))
      (if (= biogrid 1)
        (set! tmp (append tmp (pathway-gene-interactors node))))
        (if (null? tmp)
          '()
          tmp
        )
      )
      ) pw)) )


  (if (equal? prot "True")
    (set! pw (append pw (find-protein (GeneNode gene) 0))) ;; when proteins are selected, genes should only be linked to proteins not to pathways
  )

(smpdb-ctr #:enter? #f)
  (append pw ls)
))

;; From reactome

(define (reactome gene prot sm pwlst go biogrid)
(reactome-ctr #:enter? #t)
    (let (
      [pw (find-pathway-member (GeneNode gene) "R-HSA")]
      [ls '()]
      )

      (set! ls (flatten (map (lambda (path)
        (let (
            [node (cog-outgoing-atom (cog-outgoing-atom path 0) 1)]
            [tmp '()]
        )
          (set! pwlst (append pwlst (list node)))
          (set! tmp (append tmp (find-pathway-genes node go)))
          (if (equal? prot "True")
            (let ([prots (cog-outgoing-set (find-mol node "Uniprot"))])
              (if (not (null? prots))
                (set! tmp (append tmp prots))
                (set! tmp (append tmp (node-info node)))))
            )
          (if (= biogrid 1)
            (set! tmp (append tmp (pathway-gene-interactors node))))
          (if (equal? sm "True")
            (set! tmp (append tmp (cog-outgoing-set (find-mol node "ChEBI"))))
          )
          (set! tmp (append tmp (list (pathway-hierarchy node pwlst))))
          (if (null? tmp)
            '()
            tmp
          )
        )
      )
      pw)))

    (if (equal? prot "True")
    (set! pw (append pw (find-protein (GeneNode gene) 1)))
    )
(reactome-ctr #:enter? #f)
      (list (append pw ls) pwlst)
  ))
