;;; MOZI-AI Annotation Scheme
;;; Copyright © 2019 Abdulrahman Semrie
;;; Copyright © 2019 Hedra Seid
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

(define-module (annotation functions)
    #:use-module (annotation util)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog bioscience)
;    #:use-module (rnrs base)
    #:use-module (srfi srfi-1)
    #:use-module (ice-9 match)
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


(define-public find-pathway-genes-ctr (accum-time "find-pathway-genes"))
(define-public add-pathway-genes-ctr (accum-time "add-pathway-genes"))
(define-public find-go-term-ctr (accum-time "find-go-term"))
(define-public find-memberln-ctr (accum-time "find-memberln"))
(define-public add-go-info-ctr (accum-time "add-go-info"))
(define-public find-parent-ctr (accum-time "find-parent"))

(define-public match-gene-interactors-ctr (accum-time "match-gene-interactors"))
(define-public find-output-interactors-ctr (accum-time "find-output-interactors"))

(define-public find-pathway-member-ctr (accum-time "find-pathway-member"))
(define-public pathway-gene-interactors-ctr (accum-time "pathway-gene-interactors"))
(define-public generate-result-ctr (accum-time "generate-result"))
(define-public generate-interactors-ctr (accum-time "generate-interactors"))
(define-public build-interaction-ctr (accum-time "build-interaction"))
(define-public pathway-hierarchy-ctr (accum-time "pathway-hierarchy"))
(define-public check-pathway-ctr (accum-time "check-pathway"))
(define-public find-protein-ctr (accum-time "find-protein"))
(define-public find-protein-form-ctr (accum-time "find-protein-form"))
(define-public find-mol-ctr (accum-time "find-mol"))
(define-public find-pubmed-id-ctr (accum-time "find-pubmed-id"))


;;Given an atom and list of namespaces finds the parents of that atom in the specified namespaces
(define-public (find-parent a b)
	(find-parent-ctr #:enter? #t)
	(let ((rv (xfind-parent a b)))
	(find-parent-ctr #:enter? #f)
	rv))
	
(define xfind-parent

  (lambda (node namespaces)
        (let (
          [atom (cog-outgoing-atom node 1)]
          [parents '()]
        )
        (for-each (lambda (ns)
          (set! parents (append parents (cog-outgoing-set (cog-execute! (BindLink
            (TypedVariable (Variable "$a") (TypeNode 'ConceptNode))
            (AndLink
              (InheritanceLink
                atom
                (VariableNode "$a"))
              (EvaluationLink
                (PredicateNode "GO_namespace")
                (ListLink
                  (VariableNode "$a")
                  (ConceptNode ns)
                )
              )
            )
          (ExecutionOutputLink
              (GroundedSchemaNode "scm: add-go-info")
                (ListLink
                  atom
                  (VariableNode "$a")
                ))
        ))))) 
        ) namespaces
      )
    parents  
  )
))

;;Finds Go terms of a gene
(define-public (find-memberln a b)
	(find-memberln-ctr #:enter? #t)
	(let ((rv (xfind-memberln a b)))
	(find-memberln-ctr #:enter? #f)
	rv))
	
(define xfind-memberln 
  (lambda (gene namespaces)
    (let ([go-atoms '()])

      (for-each (lambda (ns)
      
        (set! go-atoms (append go-atoms (cog-outgoing-set (cog-execute! (BindLink
            (TypedVariable (Variable "$a") (TypeNode 'ConceptNode))
            (AndLink
              (MemberLink
                gene
                (VariableNode "$a"))
                (EvaluationLink
                (PredicateNode "GO_namespace")
                (ListLink
                  (VariableNode "$a")
                  (ConceptNode ns)
                )
              )
            ) 
           (ExecutionOutputLink
              (GroundedSchemaNode "scm: add-go-info")
                (ListLink
                  gene
                  (VariableNode "$a")
                ))
          )))))
      ) namespaces)
          go-atoms
          ))
)

;;Add information for GO nodes
(define-public (add-go-info a b)
	(add-go-info-ctr #:enter? #t)
	(let ((rv (xadd-go-info a b)))
	(add-go-info-ctr #:enter? #f)
	rv))
	
(define-public (xadd-go-info child-atom parent-atom)
  (if (and (or (equal? (cog-type child-atom) 'GeneNode) (equal? (cog-type child-atom) 'MoleculeNode))
      (equal? (list-ref (string-split (cog-name parent-atom) #\:) 0) "GO"))
    (ListLink  
      (MemberLink
          child-atom
          parent-atom
      )
      (go-info parent-atom)
    )
    (begin
      (if (equal? (list-ref (string-split (cog-name parent-atom) #\:) 0) "GO")
        (ListLink 
            (InheritanceLink
              child-atom
              parent-atom
          )
          (go-info parent-atom)
        )
    ))
  )
)

;;the main function to find the go terms for a gene with a specification of the parents
(define-public (find-go-term a b c)
	(find-go-term-ctr #:enter? #t)
	(let ((rv (xfind-go-term a b c)))
	(find-go-term-ctr #:enter? #f)
	rv))
	
(define-public xfind-go-term 
  (lambda (g namespaces p)
; (format #t "duuude enter find-go-term g=~A na=~A p=~A\n" g namespaces p)
      (let (
        [res (find-memberln g namespaces)]   
      )
      (define parents (flatten (let loop (
        [i p]
        [ls res]
        [acc '()]
      )
      (cond 
        [(= i 0) (append ls acc)]
        [(null? ls) acc]
        [else (cons (loop (- i 1)  (find-parent (car (cog-outgoing-set (car ls))) namespaces) (append ls acc)) (loop i (cdr ls) '()))
          ]
      )
      )))
       (append (node-info g) parents)
    )
))

;; Finds go terms for a proteins coded by the given gene
(define-public find-proteins-goterm
  (lambda (gene namespace parent)
  (let ([prot (find-protein-form gene)]
       [annotation '()])
  (if (equal? (find-memberln prot namespace) '())
    (begin
      (let ([goterms (flatten (map (lambda (ns)
        (cog-outgoing-set (cog-execute! (BindLink
        (TypedVariable (VariableNode "$g") (Type 'ConceptNode))
        (AndLink (MemberLink gene (VariableNode "$g"))
          (EvaluationLink (PredicateNode "GO_namespace") (ListLink (VariableNode "$g") (ConceptNode ns)))
        )
        (VariableNode "$g"))))
      ) namespace))])
      (set! annotation (map (lambda (go)
            (MemberLink (stv 0.0 0.0) prot go)
      )goterms))
    ))
    (set! annotation (find-go-term prot namespace parent))
  )
  (ListLink
    annotation
    (node-info prot)
    (EvaluationLink (PredicateNode "expresses") (ListLink gene prot))
  )
  )
))
;; Add details about the GO term
(define (go-info go)
  (list
      (cog-outgoing-set (find-go-name go))
      (EvaluationLink 
        (PredicateNode "GO_namespace") 
        (ListLink 
          go 
          (if (equal? (find-GO-ns go) '()) (ConceptNode "") (find-GO-ns go))))
  )
)

;; Finds parents of a GO term ( of given namespace type) 
(define find-GO-ns 
  (lambda (go)
    (cog-outgoing-set (cog-execute!
            (GetLink
                (TypedVariable (Variable "$v") (TypeNode 'ConceptNode))
                (EvaluationLink 
                    (PredicateNode "GO_namespace")
                    (ListLink 
                      go
                      (VariableNode "$v")
                    )
                 )
                )
    ))
  
  )
)

;; Finds the name of a GO term
(define find-go-name
    (lambda(go)
        (cog-execute! (BindLink
           (TypedVariable (Variable "$a") (TypeNode 'ConceptNode))
            (EvaluationLink
               (PredicateNode "GO_name")
               (ListLink
               go
               (VariableNode "$a")
              )
            )
              (EvaluationLink
               (PredicateNode "GO_name")
               (ListLink
               go
               (VariableNode "$a")
              )
            )
           )
        )
    )
)

;;finds go definition for parser function
(define find-godef
    (lambda (go)
       (cog-execute!
        (BindLink
         (VariableNode "$def")

         (EvaluationLink
          (PredicateNode "GO_definition")
          (ListLink
           go
           (VariableNode "$def")
          )
         )
          (EvaluationLink
          (PredicateNode "GO_definition")
          (ListLink
           go
           (VariableNode "$def")
          )
         )
        )       
      )
     )
)

(define-public (find-pathway-member a b)
	(find-pathway-member-ctr #:enter? #t)
	(let ((rv (xfind-pathway-member a b)))
	(find-pathway-member-ctr #:enter? #f)
	rv))
	
(define-public (xfind-pathway-member gene db)
  (cog-outgoing-set (cog-execute! (BindLink
      (TypedVariable (Variable "$a") (TypeNode 'ConceptNode))
      (AndLink
;; FIXME it would be faster to just use srfi-1 filter here... OK?
        (EvaluationLink
          (GroundedPredicateNode "scm: filter-atoms")
          (ListLink 
            (VariableNode "$a")
            (ConceptNode db)
          )
         )
        (MemberLink
          gene
          (VariableNode "$a"))
      )
      (ExecutionOutputLink
              (GroundedSchemaNode "scm: add-pathway-info")
                (ListLink
                  gene
                  (VariableNode "$a")
                ))
    )))
)

(define-public add-pathway-info 
  (lambda (gene pathway)
     (let ([res '()])
      (if  (string-contains (cog-name pathway) "R-HSA")
        (set! res (ListLink 
            (MemberLink gene pathway)
            (node-info pathway)
        ))
     )
     (if (string-contains (cog-name pathway) "SMP")
        (set! res (ListLink 
            (MemberLink gene pathway)
            (node-info pathway)
        ))
     )
     res
     )
  )
)
;; finds genes which codes the proteins in a given pathway 
(define-public find-pathway-genes
  (lambda (pathway go)
(find-pathway-genes-ctr #:enter? #t)
(let ((rv
    (cog-outgoing-set (cog-execute! (BindLink
      (VariableList 
        (TypedVariable (VariableNode "$p") (Type 'MoleculeNode))
        (TypedVariable (VariableNode "$g") (Type 'GeneNode)))
      (AndLink
      (MemberLink
        (VariableNode "$p")
        pathway)
      (EvaluationLink
        (PredicateNode "expresses")
          (ListLink
            (VariableNode "$g")
            (VariableNode "$p") )))
      (ExecutionOutputLink
      (GroundedSchemaNode "scm: add-pathway-genes")
        (ListLink
          pathway
          (VariableNode "$g")
          go
        ))
  )))
))
(find-pathway-genes-ctr #:enter? #f)
rv)
))

(define-public (add-pathway-genes pathway gene go)
	(add-pathway-genes-ctr #:enter? #t)
	(let ((rv (xadd-pathway-genes pathway gene go)))
	(add-pathway-genes-ctr #:enter? #f)
	rv))
	
(define-public (xadd-pathway-genes pathway gene go)
; (format #t "duuude in ad-pathway path=~A gene=~A go=~A\n" pathway gene go)
(if (null? (cog-outgoing-set go))
  (ListLink
        (MemberLink gene pathway)
        (node-info gene)
        (locate-node gene)
  )
  (begin
  (let ([namespace (car (cog-outgoing-set go))]
        [parent (cadr (cog-outgoing-set go))])
  (ListLink
        (MemberLink gene pathway)
        (node-info gene)
        (locate-node gene)
        (ListLink (ConceptNode "gene-go-annotation") (find-go-term gene  (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
        (ListLink (ConceptNode "gene-pathway-annotation"))
        )
  )))
))
;; Finds proteins a gene expresses
(define-public (find-protein a b)
	(find-protein-ctr #:enter? #t)
	(let ((rv (xfind-protein a b)))
	(find-protein-ctr #:enter? #f)
	rv))
	
(define-public xfind-protein
    (lambda (gene option)
        (cog-outgoing-set (cog-execute! (BindLink
          (VariableList
            (TypedVariable (Variable "$a") (TypeNode 'MoleculeNode))
            (TypedVariable (Variable "$pw") (TypeNode 'ConceptNode)))
           (AndLink
            (MemberLink
              gene
              (VariableNode "$pw"))
            (MemberLink
              (VariableNode "$a")
              (VariableNode "$pw"))
            (EvaluationLink
            (PredicateNode "expresses")
              (ListLink
                gene
                (VariableNode "$a") ))
            )
        (ExecutionOutputLink
          (GroundedSchemaNode "scm: filter-pathway")
            (ListLink
              gene
              (VariableNode "$a")
              (VariableNode "$pw")
              (Number option)
            )
        )
      )
    )
  )
))

(define-public filter-pathway (lambda (gene prot pathway option)
  (if (and (string=? (find-prefix prot) "Uniprot") )
    (cond ((and (string-contains (cog-name pathway) "SMP") (equal? option (Number "0")))
    (ListLink
      (EvaluationLink
        (PredicateNode "expresses")
          (ListLink
            gene
            prot ))
    (node-info pathway)
    ))
    ((and (equal? option (Number "1")) (string-contains (cog-name pathway) "R-HSA"))
    (ListLink
      (EvaluationLink
        (PredicateNode "expresses")
          (ListLink
            gene
            prot ))
      (node-info pathway)
      (ListLink 
        (add-loc (MemberLink gene pathway))
      )
    )))
)))

(define (find-prefix node)
  (if (equal? (length (string-split (cog-name node) #\:)) 1)
        (cog-name node)
        (car  (string-split (cog-name node) #\:))
    )
)
;; Find heirarchy of the reactome pathway
(define-public (pathway-hierarchy a b)
	(pathway-hierarchy-ctr #:enter? #t)
	(let ((rv (xpathway-hierarchy a b)))
	(pathway-hierarchy-ctr #:enter? #f)
	rv))

(define-public (run-query QUERY)
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

(define-public (xpathway-hierarchy pw lst)
	; (format #t "Enter pathway-hierarchy sizeof lst: ~A\n" (length lst))
	(if (not (member pw lst)) '()
	(let* (
			[parents (run-query
				(Get (Variable "$parentpw")
					(Inheritance pw (Variable "$parentpw"))))]
			; [junk (format #t "pathway-hierarchy found parents=~A\n" (length parents))]
			[res-parent (map
					(lambda (parent-pw) (check-pathway pw parent-pw lst))
					parents)]
			[childs (run-query
				(Get (Variable "$childpw")
					(Inheritance (Variable "$childpw") pw)))]
			; [jank (format #t "pathway-hierarchy found childs=~A\n" (length childs))]
			[res-child (map
					(lambda (child-pw) (check-pathway child-pw pw lst))
					childs)]
		)
		(append res-parent res-child)
)))


(define-public oldpathway-hierarchy
  (lambda (pw lst)
    (let ([res-parent
      (cog-outgoing-set (cog-execute! (BindLink
        (VariableNode "$parentpw")
          (InheritanceLink
            pw
          (VariableNode "$parentpw"))
        (ExecutionOutputLink
          (GroundedSchemaNode "scm: check-pathway")
          (ListLink
            pw
            (VariableNode "$parentpw")
            (ListLink lst)
          )
        ))
      ))
    ]
    [res-child (cog-outgoing-set (cog-execute! (BindLink
      (VariableNode "$parentpw")
      (InheritanceLink
        (VariableNode "$parentpw")
        pw)
      (ExecutionOutputLink
        (GroundedSchemaNode "scm: check-pathway")
        (ListLink
         (VariableNode "$parentpw")
         pw
         (ListLink lst)
        )
      ))
    ))]
  )
  (append res-parent res-child)
)))

(define-public (check-pathway a b c)
	(check-pathway-ctr #:enter? #t)
	(let ((rv (xcheck-pathway a b c)))
	(check-pathway-ctr #:enter? #f)
	rv))
	
(define-public (xcheck-pathway pw parent-pw lst)
	(if (member pw lst) (Inheritance pw parent-pw)))

(define-public oldcheck-pathway
  (lambda (pw parent-pw lst)
    (if (and (member parent-pw (cog-outgoing-set lst)) (member pw (cog-outgoing-set lst)))
    (ListLink
      (InheritanceLink
      pw
      parent-pw)
    ))
))

;; Finds molecules (proteins or chebi's) in a pathway 
(define-public (find-mol a b)
       (find-mol-ctr #:enter? #t)
	(let ((rv (xfind-mol a b)))
	(find-mol-ctr #:enter? #f)
	rv))
	
(define-public (xfind-mol path identifier)
  (cog-execute! (BindLink
    (TypedVariable (Variable "$a") (TypeNode 'MoleculeNode))
    (AndLink
      (EvaluationLink
        (GroundedPredicateNode "scm: filter-atoms")
        (ListLink
          (VariableNode "$a")
          (ConceptNode identifier)
        )
      )
       (MemberLink
       (VariableNode "$a")
       path)
    )
    (ExecutionOutputLink
      (GroundedSchemaNode "scm: add-mol-info")
      (ListLink
        (VariableNode "$a")
        path
      )
    )))
)

;; Find coding Gene for a given protein
(define-public find-coding-gene
  (lambda (protein)
  (cog-outgoing-set (cog-execute! (BindLink
    (TypedVariable (Variable "$g") (TypeNode 'GeneNode))
    (EvaluationLink
      (PredicateNode "expresses")
      (ListLink
        (VariableNode "$g")
        protein
      )
    )
    (EvaluationLink
      (PredicateNode "expresses")
      (ListLink
        (VariableNode "$g")
        protein
      )
    )
  )
  )
)))

(define-public add-mol-info
  (lambda (mol path)
  (if (string-contains (cog-name path) "R-HSA")
    (ListLink
      (MemberLink mol path)
      (if (string-contains (cog-name mol) "Uniprot")
        (find-coding-gene mol)
        '()
        )
      (node-info mol)
      (ListLink 
        (add-loc (MemberLink mol path))
      )
    )
    (ListLink
      (MemberLink mol path)
      (if (string-contains (cog-name mol) "Uniprot")
        (find-coding-gene mol)
        '()
      )
      (node-info mol)
      (ListLink (locate-node mol))
    )
  )
))

(define-public filter-atoms
  (lambda (atom identifier)
    (if (string-contains (cog-name atom) (cog-name identifier))
        (cog-new-stv 1 1)
        (cog-new-stv 0 0) 
    )
  )
) 


;; Finds genes interacting with a given gene
(define-public (match-gene-interactors a b c)
	(match-gene-interactors-ctr #:enter? #t)
	(let ((rv (xmatch-gene-interactors a b c)))
	(match-gene-interactors-ctr #:enter? #f)
	rv))
	
; FIXME search can be made easier by using SetLink for "interacts_with"
(define-public xmatch-gene-interactors
    (lambda (gene prot go)
; (format #t "duuude match-gene=~A prot=~A go=~A\n" gene prot go)
        (cog-outgoing-set (cog-execute! (BindLink
            (VariableList
            (TypedVariable (VariableNode "$a") (Type 'GeneNode)))
              (ChoiceLink 
                (EvaluationLink
                  (PredicateNode "interacts_with")
                  (ListLink
                  gene
                  (VariableNode "$a")
                  )
                )
                (EvaluationLink
                  (PredicateNode "interacts_with")
                  (ListLink
                   (VariableNode "$a")
                   gene
                  )
                )       
             )
            
            (ExecutionOutputLink
              (GroundedSchemaNode "scm: generate-result")
                (ListLink
                  gene
                  (VariableNode "$a")
                  (Number prot)
                  go
                ))        
            )
        )))	
)

;;; Finds output genes interacting eachother 
(define-public (find-output-interactors a b c)
	(find-output-interactors-ctr #:enter? #t)
	(let ((rv (xfind-output-interactors a b c)))
	(find-output-interactors-ctr #:enter? #f)
	rv))
	
(define-public xfind-output-interactors
    (lambda(gene prot go)
        (cog-outgoing-set 
          (cog-execute! (BindLink
          (VariableList
            (TypedVariable (VariableNode "$a") (Type 'GeneNode))
            (TypedVariable (VariableNode "$b") (Type 'GeneNode)))

          (And  
            (EvaluationLink
               (PredicateNode "interacts_with")
               (ListLink
               gene
               (VariableNode "$a")
              ))

            (EvaluationLink
               (PredicateNode "interacts_with")
               (ListLink
                (VariableNode "$a")
                (VariableNode "$b")
            ))

            (EvaluationLink
               (PredicateNode "interacts_with")
               (ListLink
                gene
               (VariableNode "$b")
              ))
          )
          (ExecutionOutputLink
            (GroundedSchemaNode "scm: generate-result")
              (ListLink
                (VariableNode "$a")
                (VariableNode "$b")
                (Number prot)
                go
              ))
        )))	
))

;; Gene interactors for genes in the pathway
(define-public pathway-gene-interactors
	(make-afunc-cache do-pathway-gene-interactors))

(define-public (do-pathway-gene-interactors a)
	(pathway-gene-interactors-ctr #:enter? #t)
	(let ((rv (xpathway-gene-interactors a)))
	(pathway-gene-interactors-ctr #:enter? #f)
	rv))
	
(define-public xpathway-gene-interactors 
  (lambda (pw)
  (cog-outgoing-set (cog-execute! (BindLink
    (VariableList
     (TypedVariable (VariableNode "$g1") (Type 'GeneNode))
     (TypedVariable (VariableNode "$g2") (Type 'GeneNode))
     (TypedVariable (VariableNode "$p1") (Type 'MoleculeNode))
     (TypedVariable (VariableNode "$p2") (Type 'MoleculeNode)))
   (AndLink
     (MemberLink (VariableNode "$p1") pw)
     (MemberLink (VariableNode "$p2") pw)
     (EvaluationLink (PredicateNode "expresses") (ListLink (VariableNode "$g1") (VariableNode "$p1")))
     (EvaluationLink (PredicateNode "expresses") (ListLink (VariableNode "$g2") (VariableNode "$p2")))
     (EvaluationLink (PredicateNode "interacts_with") (ListLink (VariableNode "$g1") (VariableNode "$g2")))
   )
  (ExecutionOutputLink
    (GroundedSchemaNode "scm: generate-interactors")
		  (ListLink
        pw
        (VariableNode "$g1")
		    (VariableNode "$g2")
		  ))
  ))
)))

(define-public (find-protein-form a)
	(find-protein-form-ctr #:enter? #t)
	(let ((rv (xfind-protein-form a)))
	(find-protein-form-ctr #:enter? #f)
	rv))
	
(define-public xfind-protein-form
  (lambda (gene)
  (let ([prot
  (cog-outgoing-set (cog-execute! (BindLink
    (VariableList
      (TypedVariable (VariableNode "$p") (Type 'MoleculeNode))
      (TypedVariable (VariableNode "$b") (Type 'ConceptNode)))
    (AndLink
      (EvaluationLink (PredicateNode "expresses") (ListLink gene (VariableNode "$p")))
      (EvaluationLink (PredicateNode "has_biogridID") (ListLink (VariableNode "$p") (VariableNode "$b")))
      (EvaluationLink (PredicateNode "has_biogridID") (ListLink gene (VariableNode "$b")))
    )
    (VariableNode "$p")
  )))])
  (if (not (null? prot))
    (car prot)
    (ListLink)
  )
  ))
)

;; Grounded schema node to add info about matched variable nodes

(define-public (generate-result a b c d)
	(generate-result-ctr #:enter? #t)
	(let ((rv (xgenerate-result a b c d)))
	(generate-result-ctr #:enter? #f)
	rv))

(define-public reported-genes (make-atom-set))	
(define-public reported-grid-pairs (make-atom-set))	

(define-public (xgenerate-result gene-a gene-b prot go)
      (let* (
            [output (find-pubmed-id gene-a gene-b)]
				[res (filter
					(lambda (ggg) (not (reported-genes ggg)))
					(list gene-a gene-b))]
				[pairs (reported-grid-pairs (Set gene-a gene-b))]

            [interaction (if (= 1 (string->number (cog-name prot))) 
                (ListLink
                  (build-interaction gene-a gene-b output "interacts_with")
                  (build-interaction (find-protein-form gene-a) (find-protein-form gene-b) output "inferred_interaction"))
                (build-interaction gene-a gene-b output "interacts_with"))]
            [namespace (if (null? (cog-outgoing-set go)) '() (car (cog-outgoing-set go)))]
            [parent (if (null? (cog-outgoing-set go)) '() (cadr (cog-outgoing-set go)))]
          )
          (if (null? interaction)
            (ListLink)
          )
          (match res
              ((a b)
                  (begin 
                      (if (= 1 (string->number (cog-name prot)))
                        (let ([coding-prot-a (find-protein-form a)]
                              [coding-prot-b (find-protein-form b)])
                        (if (or (equal? coding-prot-a (ListLink)) (equal? coding-prot-b (ListLink)))
                          (ListLink)
                          (ListLink
                            interaction
                            (Evaluation (Predicate "expresses") (ListLink a coding-prot-a))
                            (node-info a)
                            (node-info coding-prot-a)
                            (locate-node coding-prot-a)
                            (EvaluationLink (PredicateNode "expresses") (ListLink b coding-prot-b))
                            (node-info b)
                            (node-info coding-prot-b)
                            (locate-node coding-prot-b))  ;;; <<< bug fix here
                        ))
                      (ListLink
                          interaction
                          (node-info a)
                          (locate-node  a)
                          (node-info b)
                          (locate-node  b)
                          (if (not (null? namespace))
                          (ListLink
                            (ConceptNode "gene-go-annotation")
                            (find-go-term a (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
                            (find-go-term b (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
                            (ListLink (ConceptNode "biogrid-interaction-annotation"))
                            )
                            '() )
                      )
                  ))
              )
              ((a)
                  (begin 
                      (if (= 1 (string->number (cog-name prot)))
                        (let ([coding-prot (find-protein-form a)])
                        (if (equal? coding-prot (ListLink))
                          (ListLink)
                          (ListLink
                            interaction
                            (EvaluationLink (PredicateNode "expresses") (ListLink a coding-prot))
                            (node-info a)
                            (node-info coding-prot)
                            (locate-node coding-prot))
                        ))
                      (ListLink
                          interaction
                          (node-info a)
                          (locate-node  a)
                          (if (not (null? namespace))
                          (ListLink (ConceptNode "gene-go-annotation")
                             (find-go-term a
                             (string-split (cog-name namespace) #\ )
                             (string->number (cog-name parent)))
                          (ListLink (ConceptNode "biogrid-interaction-annotation"))
                          )
                          '()
                          ))
                  ))
              )
              (()
                  (if pairs
                    (ListLink)
                    (ListLink
                          interaction
                      )
                  )
              )
          )
      )
)
    
(define-public (oldgenerate-result gene-a gene-b prot go)

    (if  
      (and (not (equal? (cog-type gene-a) 'VariableNode)) (not (equal? (cog-type gene-b) 'VariableNode))
        )  
      (let* (
            [output (find-pubmed-id gene-a gene-b)]
            [res (flatten (map (lambda (x) 
                              (if (not (member (cog-name x) (biogrid-genes)))
                                  (cog-name x)
                                  '()
                              ) 
              )  (list gene-a gene-b))) ]
            [interaction (if (= 1 (string->number (cog-name prot))) 
                (ListLink
                  (build-interaction gene-a gene-b output "interacts_with")
                  (build-interaction (find-protein-form gene-a) (find-protein-form gene-b) output "inferred_interaction"))
                (build-interaction gene-a gene-b output "interacts_with"))]
            [namespace (if (null? (cog-outgoing-set go)) '() (car (cog-outgoing-set go)))]
            [parent (if (null? (cog-outgoing-set go)) '() (cadr (cog-outgoing-set go)))]
            [pairs (find (lambda (x) (equal? x (cons (cog-name gene-a) (cog-name gene-b)))) (biogrid-pairs))]
          )
          (if (null? interaction)
            (ListLink)
          )
          (if (not pairs)
            (biogrid-pairs (append (biogrid-pairs) (list (cons (cog-name gene-a) (cog-name gene-b)))))
          )
          (match res
              ((a b)
                  (begin 
                      (biogrid-genes (append (list a b) (biogrid-genes)))
                      (if (= 1 (string->number (cog-name prot)))
                        (let ([coding-prot-a (find-protein-form (GeneNode a))]
                              [coding-prot-b (find-protein-form (GeneNode b))])
                        (if (or (equal? coding-prot-a (ListLink)) (equal? coding-prot-b (ListLink)))
                          (ListLink)
                          (ListLink
                            interaction
                            (EvaluationLink (PredicateNode "expresses") (ListLink (GeneNode a) coding-prot-a))
                            (node-info (GeneNode a))
                            (node-info coding-prot-a)
                            (locate-node coding-prot-a)
                            (EvaluationLink (PredicateNode "expresses") (ListLink (GeneNode b) coding-prot-b))
                            (node-info (GeneNode b))
                            (node-info coding-prot-b)
                            (locate-node coding-prot-a))
                        ))
                      (ListLink
                          interaction
                          (node-info (GeneNode a))
                          (locate-node  (GeneNode a))
                          (node-info (GeneNode b))
                          (locate-node  (GeneNode b))
                          (if (not (null? namespace))
                          (ListLink
                            (ConceptNode "gene-go-annotation")
                            (find-go-term a (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
                            (find-go-term b (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
                            (ListLink (ConceptNode "biogrid-interaction-annotation"))
                            )
                            '() )
                      )
                  ))
              )
              ((a)
                  (begin 
                      (biogrid-genes (append (list a) (biogrid-genes)))
                      (if (= 1 (string->number (cog-name prot)))
                        (let ([coding-prot (find-protein-form (GeneNode a))])
                        (if (equal? coding-prot (ListLink))
                          (ListLink)
                          (ListLink
                            interaction
                            (EvaluationLink (PredicateNode "expresses") (ListLink (GeneNode a) coding-prot))
                            (node-info (GeneNode a))
                            (node-info coding-prot)
                            (locate-node coding-prot))
                        ))
                      (ListLink
                          interaction
                          (node-info (GeneNode a))
                          (locate-node  (GeneNode a))
                          (if (not (null? namespace))
                          (ListLink (ConceptNode "gene-go-annotation") (find-go-term a (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
                          (ListLink (ConceptNode "biogrid-interaction-annotation"))
                          )
                          '()
                          ))
                  ))
              )
              (()
                  (if pairs
                    (ListLink)
                    (ListLink
                          interaction
                      )
                  )
              )
          )
      )
  (ListLink)
)
    
)

(define-public (build-interaction a b c d)
	(build-interaction-ctr #:enter? #t)
	(let ((rv (xbuild-interaction a b c d)))
	(build-interaction-ctr #:enter? #f)
	rv))
	
(define-public (xbuild-interaction interactor-1 interactor-2 pubmed interaction_pred)
  (if (or (equal? (cog-type interactor-1) 'ListLink) (equal? (cog-type interactor-2) 'ListLink))
    '()
    (if (null? pubmed) 
      (EvaluationLink 
        (PredicateNode interaction_pred) 
        (ListLink interactor-1 interactor-2))
      (EvaluationLink
        (PredicateNode "has_pubmedID")
        (ListLink (EvaluationLink 
                  (PredicateNode interaction_pred) 
                  (ListLink interactor-1 interactor-2))  
                pubmed))
    )
  )
)

(define-public (generate-interactors a b c)
	(generate-interactors-ctr #:enter? #t)
	(let ((rv (xgenerate-interactors a b c)))
	(generate-interactors-ctr #:enter? #f)
	rv))
	
(define-public reported-pathways (make-atom-set))

(define-public (xgenerate-interactors path var1 var2)
	(if (reported-pathways (Set var1 var2)) '()
		(let ([output (find-pubmed-id var1 var2)])
			(if (null? output) 
				(EvaluationLink 
					(PredicateNode "interacts_with") 
					(ListLink var1 var2))
				(EvaluationLink
					(PredicateNode "has_pubmedID")
					(ListLink 
						(EvaluationLink 
							(PredicateNode "interacts_with") 
							(ListLink var1 var2))
						output)))))
)

(define-public (oldgenerate-interactors path var1 var2)
      (if (not (string=? (cog-name var1) (cog-name var2)))
      
      (let ([pairs (find (lambda (x) (or (equal? (cons (cog-name var1) (cog-name var2)) x)
                                          (equal? (cons (cog-name var2) (cog-name var1)) x)
                                      )
        
                  ) (biogrid-pairs-pathway))]
          )
          (if pairs
            '()
            (let (
              [output (find-pubmed-id var1 var2)]
              )
                (biogrid-pairs-pathway (append (biogrid-pairs-pathway) (list (cons (cog-name var1) (cog-name var2)))))
               (if (null? output) 
                (EvaluationLink 
                  (PredicateNode "interacts_with") 
                  (ListLink var1 var2))
                (EvaluationLink
                  (PredicateNode "has_pubmedID")
                  (ListLink 
                    (EvaluationLink 
                        (PredicateNode "interacts_with") 
                        (ListLink var1 var2))  
                          output)
                        )
                )
              )
            )
          )
        '()
    )
)

;;                           

(define-public (find-pubmed-id a b)
	(find-pubmed-id-ctr #:enter? #t)
	(let ((rv (xfind-pubmed-id a b)))
	(find-pubmed-id-ctr #:enter? #f)
	rv))
	
(define-public (xfind-pubmed-id gene-a gene-b)
	(cache-find-pubmed-id (Set gene-a gene-b)))
	
(define cache-find-pubmed-id
	(make-afunc-cache do-find-pubmed-id))

(define (do-find-pubmed-id gene-set)
"
	This is expecting a (SetLink (Gene \"a\") (Gene \"b\")
"
 (let* (
	[gene-a (cog-outgoing-atom gene-set 0)]
	[gene-b (cog-outgoing-atom gene-set 1)]
	[pub (run-query
     (GetLink
       (VariableNode "$pub")
       (EvaluationLink
           (PredicateNode "has_pubmedID")
           (ListLink
            (EvaluationLink 
                (PredicateNode "interacts_with") 
                  (ListLink
                    gene-a
                    gene-b
                  ))
            (VariableNode "$pub")
           )
         )

   ))])
   (if (null? pub)
     (set! pub (run-query
     (GetLink
       (VariableNode "$pub")
       (EvaluationLink
           (PredicateNode "has_pubmedID")
           (ListLink
             (EvaluationLink 
                (PredicateNode "interacts_with") 
                  (ListLink
                    gene-b
                    gene-a
                  ))
             (VariableNode "$pub")
           )
         )
   ))
   ))
   pub
))

(define-public (find-crna gene protein)
  (cog-execute! (BindLink
  (VariableList
    (TypedVariable (Variable "$a") (TypeNode 'MoleculeNode))
    (TypedVariable (Variable "$b") (TypeNode 'MoleculeNode)))
    (AndLink
      (EvaluationLink
        (PredicateNode "transcribed_to")
        (ListLink
          gene
          (VariableNode "$a")
        )
      )
      (EvaluationLink
        (PredicateNode "translated_to")
        (ListLink
          (VariableNode "$a")
          (VariableNode "$b")
        )
      )
    )
    (ExecutionOutputLink
      (GroundedSchemaNode "scm: filternc")
        (ListLink 
          gene
          (VariableNode "$a")
          (VariableNode "$b")
          (Number protein))
    )
))
)

(define-public (find-ncrna gene)
  (cog-execute! (BindLink
    (TypedVariable (VariableNode "$a") (TypeNode 'MoleculeNode))
      (EvaluationLink
        (PredicateNode "transcribed_to")
        (ListLink
          gene
          (VariableNode "$a")
        )
      )
      (ExecutionOutputLink
        (GroundedSchemaNode "scm: filternc")
		      (ListLink 
            gene
            (VariableNode "$a")
            (ListLink)
            (Number 0))
		  )
))
)

;; filter non-coding RNA
(define-public (filternc gene rna prot prot-switch)
  (if (equal? prot (ListLink))
      (ListLink
            (EvaluationLink
              (PredicateNode "transcribed_to")
              (ListLink
                gene
                rna)
            )
            (node-info rna)
      )
      (if (equal? (cog-name prot-switch) "1")
        (ListLink
          (EvaluationLink
            (PredicateNode "transcribed_to")
              (ListLink
                  gene
                  rna))
          (EvaluationLink
            (PredicateNode "translated_to")
              (ListLink
                  rna
                  prot))
          (node-info rna)
          (node-info prot)
        )
        (ListLink
          (EvaluationLink
            (PredicateNode "transcribed_to")
              (ListLink
                  gene
                  rna))
          (node-info rna)
        )
      )
  )
)
