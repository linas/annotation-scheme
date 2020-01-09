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
(define-module (annotation biogrid)
	#:use-module (annotation functions)
	#:use-module (annotation util)
	#:use-module (opencog)
	#:use-module (opencog exec)
	#:use-module (opencog bioscience)
	#:use-module (annotation parser)
	#:export (biogrid-interaction-annotation)
)
(define* (biogrid-interaction-annotation gene-nodes file-name #:key (interaction "Proteins") (namespace "") (parents 0))
  (let ([result '()]
[gctr 0]
        [go (if (string=? namespace "") (ListLink) 
                (ListLink (ConceptNode namespace) (Number parents)))])
	
; FIXME: (find-output-interactors (GeneNode gene) 1 go) is called twice,
; once for proteins, once for genes.  No need to do that.
; Also should throw error if neither is set.
	(for-each (lambda (gene)
(set! gctr (+ 1 gctr))
		(if (equal? interaction "Proteins")
(let ((start (get-internal-real-time)))
			(set! result (append result (match-gene-interactors (GeneNode gene) 1 go) (find-output-interactors (GeneNode gene) 1 go)))
(format #t "Did grid-protein ~A of ~A for ~A result-len=~A time=~6f\n"
gctr numg gene (length result) (* 1.0e-9 (- (get-internal-real-time) start)))

)
		)

		(if (equal? interaction "Genes") 
(let ((start (get-internal-real-time)))
				(set! result (append result  (match-gene-interactors (GeneNode gene) 0 go) (find-output-interactors (GeneNode gene) 0 go)))
(format #t "Did grid-gene ~A of ~A for ~A result-len=~A time=~6f\n"
gctr numg gene (length result) (* 1.0e-9 (- (get-internal-real-time) start)))
)
		)
	) gene-nodes)

	 (let (
    	[res (ListLink (ConceptNode "biogrid-interaction-annotation") (ListLink result))]
  		)
    	(write-to-file res file-name "biogrid")
		res
  	)
))
