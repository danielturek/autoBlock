__Daniel Turek__

_Stage: alpha_

~1500-2000 words (~6-9 pages)

#### BACKGROUND
_Short answers (a few sentences each) to general questions about reproducibility and scientific research_

1) Who are you and what is your research field? Include your name, affiliation, discipline, and the context of your research program necessary to frame your case study.  

My name is Daniel Turek, and my area of research is computational statistics and algorithms, with generally with application to problems of ecological statistics.  The case study I’ll be describing relates to the development of an automated procedure for improving MCMC sampling efficiency.

2) Define what the term "reproducibility" means to you generally and/or in the particular context of your case study.

In the context of my case study, “reproducibility” means that users / reviewers can re-create the same improvements in MCMC efficiency, when applying our procedure on the same example input models.  However, the results will not match exactly due to stochastic differences, but should exhibit similar order of magnitude improvements.

3) Why do you think that reproducibility in your domain is important?

Reproducibility is important so that others may verify the results stated in our publication of this procedure are genuine.

4) How or where did you learn the reproducible practices described in your case study? Mentors, classes, workshops, etc.

BIDS.  But seriously, mostly through practice of using GitHub, as well as general programming experience.  No specific classes or workshops come to mind.

5) What do you see as the major pitfalls to doing reproducible research in your domain, and do you have any suggestions for working around these? Examples could include legal, logistical, human, or technical challenges.

In the area of computational statistics, there would be few barriers to reproducible research aside from ignorance or technical inability.  However, this case study does highlight one genuine barrier: that of performance differences between various machines and computing platforms, which will affect algorithm runtime, which factors into our measure of efficiency.

6) What do you view as the major incentives for doing reproducible research?

Primarily, as stated in question 3, so that others may actually (and easily) verify our results, if they so choose.

#### WORKFLOW NARRATIVE
_Workflow narrative and diagram - a textual narrative describing your workflow with an accompanying diagram showing your workflow visually_  

The process begins with our team brainstorming how an automated procedure for improving MCMC efficiency could work.  This is arguably the most fun part of the entire process.  Anywhere from 2-4 people actually hitting the whiteboard to discuss ideas.  Each of several sessions lasts a few hours.  We review theory and literature between these sessions, too.  This initial exploration occurs over one or two weeks.

A plausible idea is hatched, and now must be prototyped to see if it actually works.  The project lead implements the procedure in R, since our engine for doing MCMC runs natively there.  This works well for our team, since everyone is comfortable in R, and code may be shared and reviewed easily.  We create a private Github repository where members of our team write/review/modify the algorithm.  This is a private repo amongst us, since it’s entirely experimental at this point, and not intended for the public.  There is little (or no) documentation at this point.

Multiple iterations are possible at this stage, whereby ideas are implemented and tested.  Depending on the results of each iteration, we go back to the drawing board several times, to figure out where the previous algorithm failed, and how it can be improved.  Once again, we implement a newer version, and test it using a small number of tests we’ve devised.  This part of the process is time consuming, and potentially frustrating, as many dead-ends are reached.  The path forward is not always clear.  This process of revising and testing our MCMC procedure may take up to 6 months.

Eventually (ideally), this process converges to a functional algorithm.  All members of our team are satisfied with the results, and agree the algorithm is ready for a more professional implementation, and hopefully publication.

One or two team members (those closest with the MCMC engine) do a more formal implementation of our procedure for improving MCMC efficiency.  This implementation is added to an existing public Github repository, which contains the basis of the MCMC engine for public use.  This step should only take a few weeks, since the algorithm is well-defined and finalized.  Appropriate documentation is also written in the form of R help files, which are also added to the public repo.

Concurrently, a suite of reproducible examples is assembled.  These come from known, standard, existing models & data sets, and also a selection of custom models, which are chosen as being either “common” applications of MCMC, or “difficult” applications of MCMC.  A new public Github repository is created, and these example models and data sets are added in the form of R data files.  Additionally, bash scripts for running the MCMC procedure on these examples are added, and a sensible README.  The sole purpose of this repo is to be referenced in the our manuscript, as a place containing fully automated scripts for reproducing the results presented in the manuscript.

The reproducible examples work well, and by including a fixed RNG seed in the executable scripts, we can guarantee the same sampling results from each MCMC.  However, the exact *timing* of each MCMC run will vary (between runs and computing platforms), and hence the final measure of efficiency will vary, too.  Thus, the *exact* results are not perfectly reproducible, but will vary approximately 5% between runs.

Team members jointly contribute to drafting a manuscript describing this new procedure, and presenting results for the suite of example models.  The manuscript specifically references the aforementioned repository, and also explains the caveat in exact reproduction of the results — namely, that they will vary slightly from those presented, and why.  The reviewers are nonetheless thrilled with the algorithm and reproducible results, and readily accept the manuscript for publication.

PAIN POINT(S)
3. Pain point(s) - an in depth discussion of ideally one, or possibly a few, sections of your workflow that you would describe as failed, incomplete, or particularly challenging in the context of reproducibility
Sections 3-5 are more detailed reflections on portions of your workflow and should each be approximately 200-400 words.

Describe in detail the steps of a reproducible workflow which you consider to be particularly painful. How do you handle these? How do you avoid them?

KEY BENEFIT(S)
4. Key achievement(s) - an in depth discussion of ideally one, or possibly a few, sections of your workflow that you are particularly proud of and feel are particularly important or useful to others

Flesh out a key step that others could benefit from.

TOOLBOX [OPTIONAL]
5. Toolbox [Optional] - if applicable, a detailed description of a particular specialized tool that plays a key role in making your workflow reproducible, if you think that the tool might be of broader interest or relevance to a general

If there is a specific tool that is critical to the achievement of reproducibility in your case study and that you feel may be underappreciated or underused, flesh it out here for benefit of others. Common tools like git or LaTeX should not be described here unless you have used them in an unconventional or particularly effective way that you feel other scientists would not know about.
