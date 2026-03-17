# Response to reviewers 

> We would like to sincerely thank the reviewers and the editor for their exceptionally fast and constructive feedback. The promptness of this review process is greatly appreciated and has allowed us to improve the manuscript significantly.

## Reviewer #2:
 The manuscript presents a geometric multigrid preconditioner for high-order SBM discretizations, introducing a novel "shy" patch smoother to effectively address the algebraic challenges induced by the inherent features of the SBM. The approach is demonstrated to yield robust and mesh-independent convergence up to polynomial degree p = 3 in both two- and three-dimensional settings. The paper tackles a timely and relevant problem, and will be interesting for the unfitted mesh and multigrid communities. The results indicate the potential of the contribution.

 However, the current presentation lacks sufficient methodological detail and clarity to meet the standards of CMAME. In particular, the description of the smoother and multigrid components is in my opinion too concise, limiting reproducibility and accessibility. I strongly encourage the authors to substantially expand and clarify the methodological section, including more detailed explanations, algorithmic descriptions, and illustrative figures, in order to make the paper more didactic and self-contained. I am convinced that with these improvements, the manuscript could become a strong and valuable contribution to the field.
 
 > *limiting reproducibility* We have released the source code.
 >
 > Added a short, self-contained subsection "Algorithm for Patch Generation" with two algorithm floats (main routine and `CollectFullResidualDoFs`) that follow the implementation.


### Main comments & suggestions
1. When referring to high-order SBM the authors may want to refer to the recent advances of Antonelli in IGA SBM. I think that the work of the authors might be of impact to this community.
    > We thank the reviewer for pointing out this relevant work. We have added a citation to the recent advances by Antonelli regarding IGA SBM in the introduction. 

2. In the introduction, I think that the sentences from "Level sets, which represent the domain" to "with corners was analyzed in [5].".
    > We have merged these duplicate sentences in the introduction to clear up the redundancy.


3. In the introduction, after having reviewed the SBM and multigrid literature, which I acknowledge is precise, the authors include two paragraphs about the CutFEM, which I think add no value to the discussion. If the authors really consider it necessary, I think it should appear before the core of the section.
    > We have compacted the comparison with CutFEM and moved it earlier in the introduction to provide a more concise context for our method.

4. When describing unfitted methods the authors state that "cells of $\mathcal{T}_h$ are classified as active based on their intersection with $\Omega$." This is customary but not always true (think on CutFEM with thin-walled bodies in which all cells are active). Based on this, I'd suggest to subtly change the statement.
    > We agree, and have changed the statement.

5. I'd suggest to refer to the $\lambda = 0.5$ as the optimal surrogate boundary approach. Besides, the authors say that this may lead to ill conditioning. I'd like to ask the authors about this statement and to be more specific about it in the manuscript.
    > We have updated the text to refer to $\lambda=0.5$ as the "optimal surrogate boundary approach". Furthermore, we have refined the statement regarding ill-conditioning. Instead of indicating formal ill-posedness, we clarified that this behavior is linked to the appearance of complex eigenvalues with non-zero imaginary parts for higher polynomial degrees associated with negative shifts, as revealed by the 1D experiments detailed in [34].

6. I would provide a more rigorous description of the high-order Taylor expansion considering that the paper targets order higher than one. Besides, I'd also include the truncation error in the definition.
    > We have refined the definition of the high-order Taylor expansion in the manuscript. The formulation now explicitly includes the polynomial summation form up to degree $k$ and formally describes the $\mathcal{O}(\|\mathbf{d}\|^{k+1})$ truncation error, ensuring a rigorous mathematical description of the boundary condition extrapolation.

7. In Eq. 3, shouldn't the penalty Nitsche terms (those with $\sigma$) be tested against the expansion of the test function? (and not the test function as in current version of the manuscript).
    > We thank the reviewer for identifying this oversight. We have corrected Equations 3 and 4; the penalty terms are now correctly tested against the extrapolation of the test function, $\mathcal{E}v_h$.

8. In Eq. 4, the notation for the penalty constant changes from $\sigma$ to $\sigma_\Gamma$. Please keep consistency. Besides, the authors define again some terms already defined before.
    > We have standardized the notation, using $\sigma$ for the penalty parameter throughout the manuscript, and eliminated redundant definitions.

9. As far as I understand from the paragraph right after Eq. 4, the values in the true boundary $\Gamma$ are computed by direct evaluation in order to skip high order derivatives. Does this mean that the values are interpolated in $\Gamma$ by using the discretization of the intersected elements? If so, this would imply to add all the DOFs of the intersected cells (not only those attached to the surrogate boundary) to the linear system but without performing the integration of such cells, something that, if I am not wrong, would result in zero entries in the left hand side matrix and thus in an ill-conditioned problem. I think that the authors should better explain this.
    > We have added a clarification to the implementation description after Equation 4. The values at the true boundary are indeed computed using the discretization of the cells within the surrogate domain. For each point on the surrogate boundary, we find its closest point on the true boundary and evaluate the shape functions at that point. Importantly, since the shape functions are polynomials defined globally per element, this evaluation is performed by mapping the physical point in the background mesh to the reference unit cell of the respective active element—meaning we are essentially computing shape function values at points that may lie slightly outside the standard $[0,1]^d$ unit cell. This approach avoids the need for explicit high-order derivative computations while maintaining the efficiency of the method, and does not require adding DOFs for excluded intersected cells to the linear system.

10. In Figure 2, I think that some of the intersected cells must be active according to the $\lambda$ criterion.
> 
<!-- We should add shading or a caption note to clarify which intersected cells are active vs inactive. -->

11. Section 3.2. What is a central vertice? Indeed, I think that this section will greatly benefit from some figures to support the discussion.
    > 
<!-- Provide a defining sentence for "central vertex" and consider adding a supporting 2D schematic diagram of a patch centered on a vertex. -->

12. Section 4.1. Which is the motivation to add all the cells' DOFs to the linear system and get an ill-conditioned problem. I would understand paying this price in the case of moving boundaries to avoid resetting the sparse matrix graph each time $\Gamma$ changes but I do not see the point in this case. Can the authors clarify this?
    > 
<!-- Clarify. We keep zero rows in the sparse matrix for code simplicity. Explain the exact mechanism.
This would be an issue in case of direct solver, and indeed we had to address it on the coarsest level where direct solver is used. This is resolved by contraining those degrees of freedom. On the fine level, we have iterative solver (GMRes) and since both the RHS and the matrix are zeroes for those rows it resutls in passing around some zoroes, but does not impact the convergence. Itreative solvers are known to handle some singular systems.

This also allowed us to just reuse the existing multigrid transfers so that we only had to implement the smoother. Since transfers can be heavily optimized reusing existing methods is important.

In an optimized implementation the number of those inactive cells could be greatly reduced by refining only cells that are at least partially inside the domain. 

 -->

13. Section 4.2. The idea described in "The extrapolation of the function values..." is not clear at all for me. I think this is very much related to my observation 9. Can the authors clarify this.
    > We have rephrased the respective paragraph in Section 4.2 to clarify the implementation details of the extrapolation mechanism, aligning it with the response to Comment 9. We now explicitly state that the extension is performed by evaluating shape functions at points projected onto the true boundary $\Gamma$, which may lie outside the standard unit cell of the active elements. This provides a direct polynomial extrapolation without the need for manual derivative computations.

14. In the numerical results, the authors set a computational domain ranging from -1.01 to 1.01 in each direction. I think this is so to avoid zero values of the level set describing the surface. This is something likely happening in industrial applications, so I would state the reasoning behind this choice as this is not resolved in the author's paper.
    > The background mesh covering $[-1.01, 1.01]$ ensures that the closest point projection onto the true boundary $\Gamma$ can be reliably computed for all points on the surrogate boundary $\tilde{\Gamma}$ within our level set framework. We have clarified this in the text.
<!-- @Copilot: Resolve that and write a short response after `>` above: The 1.01 was chosen to ensure that the true boundary lies within the mesh: this is because we are using level set. -->
15. Figure 3 (right). What do the dashed lines represent?
    > The dashed lines in Figure 3 (right) indicate the maximum negative shift magnitude, multiplied by $-1$ for consistency in the visualization. 
<!-- Add a trace in the caption specifying that dashed lines denote theoretical/expected slopes or references, if applicable. -->

16. Though explained in [34], I think that the readers would appreciate a more detailed remark on the fact that using the optimal surrogate boundary (i.e., $\lambda=0.5$) results in negative eigenvalues that affect the convergence as, a priori, one may easily think that $\lambda=0.5$ is always the best.
    > We have expanded our discussion on this phenomenon in the numerical results section. As you pointed out, while it is intuitive to expect the optimal surrogate boundary geometrically ($\lambda=0.5$) to perform best, doing so leads to shifts that point inwards to the true domain. Following our earlier investigation (Ref [34]), we explicitly clarified in the text that these negative shifts result in complex eigenvalues with non-zero *imaginary* (not negative) parts for higher-order elements, which disrupts the smoothing properties of the multigrid preconditioner and causes failures.

### Minor comments
- Use Latex format quotes (i.e., `` '')
    > We have updated the quotes to emphasize phrasing as suggested.  
<!-- @Copilot: resolve that, reply to the reviewer. Also in general prefer \emph instead of qoutes. -->
- The authors define the acronyms several times. This happens recursively for Continuous Galerkin, Discontinuous Galerkin, Shifted Boundary Method, Degrees of Freedom, etc.
    > We have reviewed the manuscript and removed redundant acronym definitions.
<!-- @Copilot: Resolve that. reply to the reviewer   -->
- I think the authors missed section 4 in the paper outline at the end of the introduction.
    > Thank you for catching this omission. We have added the implementation details section to the outline.
<!-- @Copilot: resolve that, reply to the reviewer -->
- When defining the FE discretization, the authors state "triangulation consisting of quads. and hexas.". I think that discretization is more suitable than triangulation here.
    > We agree, and have updated the terminology from triangulation to discretization.

- "If the threshold is not met within 100 iterations, the solver is considered to have failed" statement is repeated.
    > We have removed the duplicate statements to avoid repetition.

- In the conclusion, I think that the authors wanted to use italics in *shyness*.
    > Indeed 
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
<!-- Add a paragraph in the methodology section drawing a brief parallel to overlapping Schwarz theory and discussing how the patch design captures high-frequency error components near the complicated surrogate boundary. -->

2. Several parameters (penalty \sigma, shyness threshold \xi, number of smoothing steps s, and cell threshold \lambda) are explored empirically. While this is reasonable, a more explicit discussion or practical guidelines for choosing these parameters in general problems would be useful for practitioners.
> 
<!-- Add a subsection or a concluding table/paragraph offering practical heuristics, e.g., default choices for these parameters based on polynomial degree p. -->

### Minor remarks
1. Page 5, correct "singifficant challenges"
    > Done.
2. Page 2, introduction, there are two almost repeated sentences:
   - Level sets,which represent the domain boundary as the zero level set of a function…
   - Level set methods, which represent the domain boundary as the zero level set of a function…
    > We have consolidated these two sentences to eliminate the repetition.
<!-- @Copilot resolve that and respond (shortly) to reviewer -->
3. Page 5, "our smoother operateS …"
    > Done.

4. In the acknowledgments: "The authors declare ..."
    > Acknowledgments updated.
