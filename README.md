# GL-AR
This is a repository for the development of GL-AR, a novel graph learning framework utilizing autoregressive coefficients to recover relationships exhibiting significant propagation delay.

## Paper Abstract
*This research presents a novel approach to graph learning, GL-AR, where estimated autoregressive coefficients play a pivotal role in facilitating the recovery of graph structures that exhibit propagation delay between connected vertices.  Graph learning involves inferring graph structures using a set of graph signals. Traditionally, graph signal processing can be used to achieve this by leveraging time-series graph signals, however the current state-of-the-art techniques utilizing this approach require relationships to lack propagation delay. This work proposes a departure from convention by employing the estimation of autoregressive coefficients of the time-series graph signals. This departure enables GL-AR to discern graph structures even in scenarios where the relationships contain significant propagation delay.*


*The work comprehensively outlines the conceptualization and implementation of GL-AR, showcasing its efficacy through applications to synthetic and real-world datasets. The results prominently highlight GL-AR's proficiency in learning graph structures from time-series graph signals, particularly when confronted with relationships characterized by nonzero propagation delays. Moreover, these results are obtained using an automated selection algorithm for problem parameters, which removes the need for the computationally intensive trial-and-error methods utilized in current state-of-the-art graph learning frameworks.*

*While underscoring the evident advantages of GL-AR, this work also conscientiously addresses inherent limitations and delves into potential directions for future research.
This multi-faceted exploration not only emphasizes the potential impact of GL-AR but additionally highlights the dynamic nature of the field and the ongoing pursuit of more robust and versatile graph learning methodologies.*

## Repository Description
This repository hosts the MATLAB codebase for GL-AR. This codebase was used to both develop and evaluate GL-AR against its predecessor, GL-SigRep. The MATLAB folder's README file goes into more detail.