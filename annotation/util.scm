;;; MOZI-AI Annotation Scheme
;;; Copyright © 2019 Abdulrahman Semrie
;;; Copyright © 2019 Hedra Seid
;;; Copyright © 2019 Enkusellasie Wondesen
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

(define-module (annotation util)
	#:use-module (opencog)
	#:use-module (opencog exec)
	#:use-module (opencog bioscience)
	#:use-module (json)
	#:use-module (ice-9 optargs)
;	#:use-module (rnrs base)
	#:use-module (rnrs exceptions)
	#:use-module (ice-9 textual-ports)
	#:use-module (ice-9 regex)
	#:use-module (srfi srfi-1)
	#:use-module (ice-9 match)
	#:export (create-node
	          create-edge filter-loc locate-node)
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


(define-public find-name-ctr (accum-time "find-name"))
(define-public locate-node-ctr (accum-time "locate-node"))
(define-public add-loc-ctr (accum-time "add-loc"))

;;Define the parameters needed for parsing and GGI
(define-public nodes (make-parameter '()))
(define-public edges (make-parameter '()))
(define-public atoms (make-parameter '()))
(define-public biogrid-genes (make-parameter '()))
(define-public biogrid-pairs (make-parameter '()))
(define-public biogrid-pairs-pathway (make-parameter '()))
(define-public annotation (make-parameter ""))
(define-public prev-annotation (make-parameter ""))

(define (get-name atom)
 (if (> (length atom) 0)
  (cog-name (car  atom))
  ""
 )
)


(define* (create-node id name defn location annotation #:optional (subgroup ""))
    (make-node (make-node-info id name defn location subgroup annotation) "nodes")
)

(define* (create-edge node1 node2 name annotation #:optional (pubmedId "") (subgroup ""))
   (make-edge (make-edge-info node2 node1 name pubmedId subgroup annotation) "edges")
)

;; Find node name and description

(define-public node-info
	(make-afunc-cache xnode-info))

(define-public (xnode-info node)
; (format #t "duude node-info=~A\n" node)
  (if (cog-node? node)
    (list
      (EvaluationLink (PredicateNode "has_name") (ListLink node (node-name node)))
    )
    (list (ListLink))
  )
)

(define (node-name node)
	(let
		( [lst (find-pathway-name node)])
; (format #t "duude in nonde-name node=~A lst=~A\n" node lst)
			(if (null? lst)
				(ConceptNode "N/A")
				(car lst)
			)
	)
)

;;Finds a name of any node (Except GO which has different structure)
(define find-pathway-name
    (lambda(pw)
			(let
        ([name 
        (if (or (string-contains (cog-name pw) "Uniprot:") (string-prefix? "ENST" (cog-name pw)))
          (let ([predicate (if (string-prefix? "ENST" (cog-name pw)) "transcribed_to" "expresses")])
            (cog-outgoing-set (cog-execute! (GetLink
              (VariableList
              (TypedVariable (VariableNode "$a") (Type 'GeneNode)))
              (EvaluationLink
                (PredicateNode predicate)
                (ListLink
                  (VariableNode "$a")
                  pw
                )
              )
            )))
          )
          (cog-outgoing-set (cog-execute! (GetLink
            (VariableList
            (TypedVariable (VariableNode "$a") (Type 'ConceptNode)))
            (EvaluationLink
                (PredicateNode "has_name")
                (ListLink
                  pw
                  (VariableNode "$a")
                )
              )	
          )))
        )])
    name
	))
)

(define-public (is-cellular-component? info-list)
	(any
		(lambda (info)
			(and
				(equal? (cog-name (gar info)) "GO_namespace")
				(equal? "cellular_component" (cog-name (gddr info)))))
		info-list)
)

(define-public (build-desc-url node)
    (cond 
        ((string-prefix? "ChEBI" node) (string-append "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=" (cadr (string-split node #\:))))
        ((string-prefix? "Uniprot" node) (string-append "https://www.uniprot.org/uniprot/" (cadr (string-split node #\:))))
        ((string-prefix? "GO" node) (string-append "http://amigo.geneontology.org/amigo/term/" node))
        ((string-prefix? "SMP" node) (string-append "http://smpdb.ca/view/" node))
        ((string-prefix? "R-HSA" node)
          (string-append "http://www.reactome.org/content/detail/" node)
        )
        ((string-prefix? "ENST" node) (string-append "http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=" node))
        ((or (string-prefix? "NM_" node) (string-prefix? "NR_" node) (string-prefix? "NP_" node) (string-prefix? "YP_" node))
          (string-append "https://www.ncbi.nlm.nih.gov/nuccore/" node))
        (else (string-append "https://www.ncbi.nlm.nih.gov/gene/"  (find-entrez (GeneNode node))))
    )
)

;; Finds entrez_id of a gene
(define (find-entrez gene)
  (let ((entrez '()))
    (set! entrez (get-name
   (cog-outgoing-set
    (cog-execute!
     (GetLink
       (VariableNode "$a")
       (EvaluationLink
        (PredicateNode "has_entrez_id")
        (ListLink
         gene
         (VariableNode "$a")
        )
       )
    )
   )
  )
  ))
   (if (equal? (length (string-split entrez #\:)) 1)
       entrez
       (cadr  (string-split entrez #\:))
   )
  )
)

(define (run-query QUERY)
"
  Call (cog-execute! QUERY), return results, delete the SetLink.
  This avoids a memory leak of SetLinks
"
   ; Run the query
   (define set-link (cog-execute! QUERY))
   ; Get the query results
   (define results (cog-outgoing-set set-link))
   ; Delete the SetLink
   (cog-delete set-link)
   ; Return the results.
   results
)

(define-public find-name
	(make-afunc-cache probe-find-name))

;;finds go name for parser function
(define (probe-find-name a)
	(find-name-ctr #:enter? #t)
	(let ((rv (oldfind-name a)))
	(find-name-ctr #:enter? #f)
	rv))
	
(define oldfind-name (lambda (atom)
     (let*
        (
          [predicate (if (regexp-match? (string-match "GO:[0-9]+" (cog-name atom))) "GO_name" "has_name")]
        )
      (get-name
       (run-query
        ; cog-execute! ; cog-execat!  ; cog-execute!
         (GetLink
          (VariableNode "$name")

          (EvaluationLink
           (PredicateNode predicate)
           (ListLink
            atom
            (VariableNode "$name")
           ))))))))


;;Given an atom and list of namespaces finds the parents of that atom in the specified namespaces
(define (newfind-name GO-ATOM)
	(define pred (Predicate
		(if (regexp-match? (string-match "GO:[0-9]+" (cog-name GO-ATOM)))
			"GO_name" "has_name")))

	(define namli (find
		(lambda (lili)
			(any
				(lambda (evli) (equal? pred (gar evli)))
				(cog-incoming-by-type lili 'EvaluationLink)))
		(cog-incoming-by-type GO-ATOM 'ListLink)))

	(if namli (cog-name (gdr namli)) "")
)

(define-public (filter-genes input-gene gene-name)
  (if (regexp-match? (string-match (string-append (cog-name input-gene) ".+$") (cog-name gene-name)))
      (stv 1 1)
      (stv 0 0) 
  )
)

(define-public (find-similar-gene gene-name)
    (cog-outgoing-set 
      (cog-execute! (BindLink
        (TypedVariable (Variable "%x") (Type "GeneNode"))
        (EvaluationLink
          (GroundedPredicateNode "scm: filter-genes")
          (ListLink
              (Gene gene-name)
              (VariableNode "%x")      
          )
        )
        (VariableNode "%x")
      ))
    
    )
)

(define (find-current-symbol gene-list)
  (let ((current (map (lambda (gene)  
    (let ((cur (cog-outgoing-set 
      (cog-execute! (BindLink
        (TypedVariable (Variable "$g") (Type "GeneNode"))
        (EvaluationLink
          (PredicateNode "has_current_symbol")
          (ListLink
              (GeneNode gene)
              (VariableNode "$g")      
          )
        )
        (VariableNode "$g")
      ))
    )))
    (if (equal? cur '()) '() (string-append gene "=>" (cog-name (car cur))))
    )
  )gene-list)))
  (flatten current)  
))

(define-public (check-outdate-genes gene-list) 
  (let ([symbol (find-current-symbol gene-list)])
    (if (equal? symbol '())
      "0"
      (string-append "1:" "The following gene symbols are outdated, here are the current " (string-join symbol ","))
    )
  )
)

(define-public (build-pubmed-url nodename)
 (string-append "https://www.ncbi.nlm.nih.gov/pubmed/?term=" (cadr (string-split nodename #\:)))
)

(define-public (write-to-file result id name)
  (catch #t (lambda ()
    (let*
        (
          [env-path (getenv "RESULT_DIR")]
          [path (if (not env-path) (string-append "/tmp/result/" id) (string-append env-path "/" id))]
          [file-name (string-append path "/" name ".scm")]
        )
        (if (not (file-exists? path))
            (mkdir path)
        )
        (call-with-output-file file-name
            (lambda (p)
              (write result p)
            )
          )
  ))  
  (lambda (key . parameters)
      (format (current-error-port) "Cannot write scheme result files. ~a: ~a\n" key parameters)
      #f
    ))
)

(define-public locate-node
	(make-afunc-cache probe-locate-node))

(define (probe-locate-node a)
   (locate-node-ctr #:enter? #t)
   (let ((rv (xlocate-node a)))
   (locate-node-ctr #:enter? #f)
   rv))

(define xlocate-node
  (lambda(node)
; (format #t "duude in locate node=~A\n" node)
      (let ([loc (run-query
        (BindLink
        (VariableNode "$go")
        (AndLink
          (MemberLink 
            node
            (VariableNode "$go"))
          (EvaluationLink
            (PredicateNode "GO_namespace")
            (ListLink
              (VariableNode "$go")
              (ConceptNode "cellular_component")))
        )
        (ExecutionOutputLink
        (GroundedSchemaNode "scm: filter-loc")
          (ListLink
            node
            (VariableNode "$go")
          ))
          )
      )
      ])
; (format #t "duude in locate loc=~A\n" loc)
      (if (null? loc)
      (set! loc 
      (run-query
        (BindLink
          (VariableNode "$loc")
          (EvaluationLink
              (PredicateNode "has_location")
              (ListLink
                node
                (VariableNode "$loc")))
          (EvaluationLink
              (PredicateNode "has_location")
              (ListLink
                node
                (VariableNode "$loc")))
          )
        ))
      )
      loc
      )
  )
)

;; filter only Cell membrane and compartments

(define-public (filter-loc node go)
; (format #t "duuude in filterloc node=~A go=~A\n" node go)
  (let ([loc (string-downcase (find-name go))])
; (format #t "duuude in filterloc loc=~A\n" loc)
  (if (or (and (not (string-contains loc "complex")) 
      (or (string-suffix? "ome" loc) (string-suffix? "ome membrane" loc))) (is-compartment loc))
        (EvaluationLink
          (PredicateNode "has_location")
          (ListLink
            node
            (ConceptNode loc)
          )
        )
  )
  ))

(define (is-compartment loc)
	(any
		(lambda (comp) (string-contains loc comp))
		(list "vesicle" "photoreceptor" "plasma" "centriole"
			"cytoplasm" "endosome" "golgi" "vacuole" "granule"
			"endoplasmic" "mitochondri" "cytosol" "peroxisome"
			"ribosomes" "lysosome" "nucle"))
)

;; Add location of a gene/Molecule node in context of Reactome pathway

(define-public (add-loc a)
   (add-loc-ctr #:enter? #t)
   (let ((rv (xadd-loc a)))
   (add-loc-ctr #:enter? #f)
   rv))

(define-public (xadd-loc node)
  (let ([child (cog-outgoing-atom node 0)] )
      (cog-outgoing-set (cog-execute!
        (BindLink
          (VariableNode "$loc")
          (AndLink
            (EvaluationLink
              (PredicateNode "has_location")
              (ListLink
                child
                (VariableNode "$loc")))
          )
            (EvaluationLink
              (PredicateNode "has_location")
              (ListLink
                child
                (VariableNode "$loc")))
          )
        )
      )
    )
)

(define-public (find-subgroup name) 
    (let ((initial (string-split name #\:)))
        (match initial
            ((a b) a )
            ((a)
                (cond 
                    ((string-prefix? "R-HSA" a) "Reactome")
                    ((string-prefix? "SMP" a) "SMPDB")
                    (else "Genes")
                )
            )
        )
    )

)
;;a helper function to flatten a list, i.e convert a list of lists into a single list
(define-public (flatten x)
  (cond ((null? x) '())
        ((and (cog-link? x) (null? (cog-outgoing-set x))) '())
        ((pair? x) (append (flatten (car x)) (flatten (cdr x))))
        (else (list x))))
