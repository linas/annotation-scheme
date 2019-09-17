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
    #:use-module (opencog query)
    #:use-module (opencog exec)
    #:use-module (opencog bioscience)
    #:use-module (rnrs base)
    #:use-module (srfi srfi-1)
    #:use-module (ice-9 match)
    
)

;;Given an atom and list of namespaces finds the parents of that atom in the specified namespaces
(define find-parent

  (lambda (node namespaces)
        (let (
          [atom (cog-outgoing-atom node 1)]
          [parents '()]
        )
        (for-each (lambda (ns)
          (set! parents (append parents (cog-outgoing-set (cog-execute! (BindLink
            (VariableNode "$a")
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
(define find-memberln 
  (lambda (gene namespaces)
    (let ([go-atoms '()])

      (for-each (lambda (ns)
      
        (set! go-atoms (append go-atoms (cog-outgoing-set (cog-execute! (BindLink
            (VariableNode "$a")
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
(define-public (add-go-info child-atom parent-atom)
  (if (and (equal? (cog-type child-atom) 'GeneNode)
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
(define-public find-go-term 
  (lambda (g namespaces p)
      (let (
        [res (find-memberln (GeneNode g) namespaces)]   
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
       (delete-duplicates parents)
       (append (node-info (GeneNode g)) parents)
    )
))

;; Add details about the GO term
(define (go-info go)
  (list
      (cog-outgoing-set (find-go-name go))
      (cog-outgoing-set (find-godef go))
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
                (VariableNode "$v")
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
            (VariableNode "$a")
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

(define-public (find-pathway-member gene db)
  (cog-outgoing-set (cog-execute! (BindLink
      (VariableNode "$a")
      (AndLink
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
(define-public (add-pathway-genes pathway gene go)
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
        (ListLink (ConceptNode "gene-go-annotation") (find-go-term (cog-name gene)  (string-split (cog-name namespace) #\ ) (string->number (cog-name parent)))
        (ListLink (ConceptNode "gene-pathway-annotation"))
        )
  )))
))
;; Finds proteins a gene expresses
(define-public find-protein
    (lambda (gene option)
        (cog-outgoing-set (cog-execute! (BindLink
          (VariableList
            (VariableNode "$a")
            (VariableNode "$pw"))
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
        )))))

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
(define-public pathway-hierarchy
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

(define-public check-pathway
  (lambda (pw parent-pw lst)
    (if (and (member parent-pw (cog-outgoing-set lst)) (member pw (cog-outgoing-set lst)))
    (ListLink
      (InheritanceLink
      pw
      parent-pw)
    ))
))

;; Finds molecules (proteins or chebi's) in a pathway 
(define-public (find-mol path identifier)
  (cog-execute! (BindLink
    (VariableNode "$a")
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
    (VariableNode "$g")
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
(define-public match-gene-interactors
    (lambda (gene prot go)
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
(define-public find-output-interactors
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
  (lambda (pw gene)
  (cog-outgoing-set (cog-execute! (BindLink
    (VariableList
      (TypedVariable (VariableNode "$g1") (Type 'GeneNode))
      (TypedVariable (VariableNode "$g2") (Type 'GeneNode))
	    (TypedVariable (VariableNode "$p") (Type 'MoleculeNode)))
    (AndLink
      (MemberLink (VariableNode "$p") pw)
      (EvaluationLink (PredicateNode "expresses") (ListLink (VariableNode "$g1") (VariableNode "$p")))
      (EvaluationLink (PredicateNode "interacts_with") (ListLink (VariableNode "$g1") (VariableNode "$g2"))) 
      (EvaluationLink (PredicateNode "interacts_with") (ListLink (VariableNode "$g2") gene))
    )
  (ExecutionOutputLink
    (GroundedSchemaNode "scm: generate-interactors")
		  (ListLink
      pw
      gene
        (VariableNode "$g1")
		    (VariableNode "$g2")
		  ))
  ))
)))

(define-public find-protein-form
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

(define-public (generate-result gene-a gene-b prot go)
    (call/cc (lambda (k)
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
                  [interaction (if (= 1 (string->number (cog-name prot))) (build-interaction (find-protein-form gene-a) (find-protein-form gene-b) output)
                                      (build-interaction gene-a gene-b output))]
                  [namespace (if (null? (cog-outgoing-set go)) '() (car (cog-outgoing-set go)))]
                  [parent (if (null? (cog-outgoing-set go)) '() (cadr (cog-outgoing-set go)))]
                  [pairs (find (lambda (x) (equal? (cons (cog-name gene-a) (cog-name gene-b)))) (biogrid-pairs))]
                )
                (if (null? interaction)
                  (k (ListLink))
                )
                (if (not pairs)
                  (biogrid-pairs (append (biogrid-pairs) (cons (cog-name gene-a) (cog-name gene-b))))
                )
                (match res
                    ((a b)
                        (begin 
                            (biogrid-genes (append (list a b) (biogrid-genes)))
                            (if (= 1 (string->number (cog-name prot)))
                              (ListLink
                                  interaction
                                  (node-info (find-protein-form (GeneNode a)))
                                  (locate-node  (find-protein-form (GeneNode a)))
                                  (node-info (find-protein-form (GeneNode b)))
                                  (locate-node  (find-protein-form (GeneNode b)))
                              )
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
                            (ListLink
                                interaction
                                (node-info (find-protein-form (GeneNode a)))
                                (locate-node  (find-protein-form (GeneNode a))))
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
    
))

(define-public (build-interaction interactor-1 interactor-2 pubmed)
  (if (or (equal? (cog-type interactor-1) 'ListLink) (equal? (cog-type interactor-2) 'ListLink))
    '()
    (if (null? pubmed) 
      (EvaluationLink 
        (PredicateNode "interacts_with") 
        (ListLink interactor-1 interactor-2))
      (EvaluationLink
        (PredicateNode "has_pubmedID")
        (ListLink (EvaluationLink 
                  (PredicateNode "interacts_with") 
                  (ListLink interactor-1 interactor-2))  
                pubmed))
    )
  )
)

(define-public (generate-ppi-result gene-a prot-a )
        (ListLink
                (EvaluationLink (PredicateNode "expresses") (ListLink gene-a prot-a))
                (node-info prot-a)
        )
)

(define-public (generate-interactors path gene var1 var2)
  (if (and (not (string=? (cog-name var1) (cog-name var2)))
          (not (or (string=? (cog-name gene) (cog-name var1))(string=? (cog-name gene) (cog-name var2))))
      )
      (let ([pairs (find (lambda (x) (equal? (cons (cog-name var1) (cog-name var2)))) (biogrid-pairs))]
          )
          (if pairs
            (ListLink)
            (let (
              [output (find-pubmed-id var1 var2)]
              )
                (biogrid-pairs (append (biogrid-pairs) (cons (cog-name var1) (cog-name var2))))
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
        (ListLink)
    )
    
)

;;                           
(define (find-pubmed-id gene-a gene-b)
 (let ([pub (cog-outgoing-set (cog-execute!
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

   )))])
   (if (null? pub)
     (set! pub (cog-outgoing-set (cog-execute!
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
   )))
   ))
   pub
))
