# Response to reviewers 
<!-- Thank reviewers and editor for insanely fast review speed -->

## Reviewer #2:
 The manuscript presents a geometric multigrid preconditioner for high-order SBM discretizations, introducing a novel "shy" patch smoother to effectively address the algebraic challenges induced by the inherent features of the SBM. The approach is demonstrated to yield robust and mesh-independent convergence up to polynomial degree p = 3 in both two- and three-dimensional settings. The paper tackles a timely and relevant problem, and will be interesting for the unfitted mesh and multigrid communities. The results indicate the potential of the contribution.

 However, the current presentation lacks sufficient methodological detail and clarity to meet the standards of CMAME. In particular, the description of the smoother and multigrid components is in my opinion too concise, limiting reproducibility and accessibility. I strongly encourage the authors to substantially expand and clarify the methodological section, including more detailed explanations, algorithmic descriptions, and illustrative figures, in order to make the paper more didactic and self-contained. I am convinced that with these improvements, the manuscript could become a strong and valuable contribution to the field.


### Main comments & suggestions
1. When referring to high-order SBM the authors may want to refer to the recent advances of Antonelli in IGA SBM. I think that the work of the authors might be of impact to this community.
> 
<!-- \cite{ANTONELLI2024shiftedIGA} -->

2. In the introduction, I think that the sentences from "Level sets, which represent the domain" to "with corners was analyzed in [5].".
> 

3. In the introduction, after having reviewed the SBM and multigrid literature, which I acknowledge is precise, the authors include two paragraphs about the CutFEM, which I think add no value to the discussion. If the authors really consider it necessary, I think it should appear before the core of the section.
> 
<!-- Shorten -->
4. When describing unfitted methods the authors state that "cells of $\mathcal{T}_h$ are classified as active based on their intersection with $\Omega$." This is customary but not always true (think on CutFEM with thin-walled bodies in which all cells are active). Based on this, I'd suggest to subtly change the statement.
> 

5. I'd suggest to refer to the $\lambda = 0.5$ as the optimal surrogate boundary approach. Besides, the authors say that this may lead to ill conditioning. I'd like to ask the authors about this statement and to be more specific about it in the manuscript.
> 

6. I would provide a more rigorous description of the high-order Taylor expansion considering that the paper targets order higher than one. Besides, I'd also include the truncation error in the definition.
> 

7. In Eq. 3, shouldn't the penalty Nitsche terms (those with $\sigma$) be tested against the expansion of the test function? (and not the test function as in current version of the manuscript).
> 

8. In Eq. 4, the notation for the penalty constant changes from $\sigma$ to $\sigma_\Gamma$. Please keep consistency. Besides, the authors define again some terms already defined before.
>   

9. As far as I understand from the paragraph right after Eq. 4, the values in the true boundary $\Gamma$ are computed by direct evaluation in order to skip high order derivatives. Does this mean that the values are interpolated in $\Gamma$ by using the discretization of the intersected elements? If so, this would imply to add all the DOFs of the intersected cells (not only those attached to the surrogate boundary) to the linear system but without performing the integration of such cells, something that, if I am not wrong, would result in zero entries in the left hand side matrix and thus in an ill-conditioned problem. I think that the authors should better explain this.
> 

10. In Figure 2, I think that some of the intersected cells must be active according to the $\lambda$ criterion.
> 

11. Section 3.2. What is a central vertice? Indeed, I think that this section will greatly benefit from some figures to support the discussion.
> 

12. Section 4.1. Which is the motivation to add all the cells' DOFs to the linear system and get an ill-conditioned problem. I would understand paying this price in the case of moving boundaries to avoid resetting the sparse matrix graph each time $\Gamma$ changes but I do not see the point in this case. Can the authors clarify this?
> 

13. Section 4.2. The idea described in "The extrapolation of the function values..." is not clear at all for me. I think this is very much related to my observation 9. Can the authors clarify this.
> 

14. In the numerical results, the authors set a computational domain ranging from -1.01 to 1.01 in each direction. I think this is so to avoid zero values of the level set describing the surface. This is something likely happening in industrial applications, so I would state the reasoning behind this choice as this is not resolved in the author's paper.
> 

15. Figure 3 (right). What do the dashed lines represent?
> 

16. Though explained in [34], I think that the readers would appreciate a more detailed remark on the fact that using the optimal surrogate boundary (i.e., $\lambda=0.5$) results in negative eigenvalues that affect the convergence as, a priori, one may easily think that $\lambda=0.5$ is always the best.
> 

### Minor comments
- Use Latex format quotes (i.e., `` '')
> 
- The authors define the acronyms several times. This happens recursively for Continuous Galerkin, Discontinuous Galerkin, Shifted Boundary Method, Degrees of Freedom, etc.
> 
- I think the authors missed section 4 in the paper outline at the end of the introduction.
> 
- When defining the FE discretization, the authors state "triangulation consisting of quads. and hexas.". I think that discretization is more suitable than triangulation here.
> 
- "If the threshold is not met within 100 iterations, the solver is considered to have failed" statement is repeated.
> 
- In the conclusion, I think that the authors wanted to use italics in *shyness*.
> 
- I think that the authors meant to write an AI declaration statement (or something of this sort) rather than Acknowledgements.
> 
- The formatting of some of the references is wrong.
> 

---

## Reviewer #3

The manuscript presents a geometric multigrid preconditioner tailored to linear systems arising from high-order Continuous Galerkin discretizations of the Shifted Boundary Method (SBM). The key contribution is the introduction of a novel Full-Residual Shy Patch smoother, a patch-based subspace correction strategy designed to overcome the well-known difficulties of standard smoothers and AMG in the presence of SBM's non-symmetric and non-local boundary coupling.

The paper is technically solid, well-motivated, and addresses an important open problem in the numerical solution of unfitted finite element methods. The proposed smoother is carefully designed, clearly explained, and supported by extensive numerical experiments in both 2D and 3D. The results convincingly demonstrate h-robustness and significantly improved performance compared to existing approaches, including AMG and previously proposed DG-SBM multigrid methods.

On the basis of that, I recommend publication after a minor revision according to the comments below.

### Major remarks
1. While the numerical evidence is strong, the paper would benefit from additional theoretical insight into why the Full-Residual Shy Patch smoother is effective. Even a qualitative spectral or error-propagation discussion (e.g., local vs. global error components, relation to Schwarz theory) would strengthen the contribution.
> 

2. Several parameters (penalty \sigma, shyness threshold \xi, number of smoothing steps s, and cell threshold \lambda) are explored empirically. While this is reasonable, a more explicit discussion or practical guidelines for choosing these parameters in general problems would be useful for practitioners.
> 

### Minor remarks
1. Page 5, correct "singifficant challenges"
> 
2. Page 2, introduction, there are two almost repeated sentences:
   - Level sets,which represent the domain boundary as the zero level set of a function…
   - Level set methods, which represent the domain boundary as the zero level set of a function…
> 
3. Page 5, "our smoother operateS …"
> 
4. In the acknowledgments: "The authors declare ..."
>
