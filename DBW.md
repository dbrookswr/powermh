---

title: 'powermh: Power analysis for studies with multiple hypotheses'

tags:
  
   - Power
   - Multiple hypotheses
   - R 

authors:

   - name: Daniel B. Wright
   
     orcid: 0000-0002-4619-1848
   
     affiliation: 1 

   - name: Marianna E. Carlucci 
     orcid: 0000-0000-0000-0000
   
     affiliation: 2
affiliations:
   - name: Alder Graduate School of Education
     index: 1
   - name: Loyola University Maryland
     index: 2 

date: 5 July 2018

bilbliography: paper.bib

---


# Summary


Power analysis is recommended by scientific societies and journals [@WilkinsonEA1999; Wright2003], and required by some funding bodies (e.g., ies.ed.gov/funding/resources.asp). Tables [@Cohen1992] and computer packages [@FaulEA2007] can assist researchers in conducting traditional power analyses. The algorithms used for traditional power analysis require the researcher to choose a single test statistic.  

The reality of many research projects is more complex. Most processes in nature cause multiple effects, and as such a useful theory of this process should make multiple predictions. Fisher summed this up nicely: ``Make your theories elaborate'' [@Cochran1965, p. 252]. Cochran explains this statement: ``What Sir Ronald meant ... was that when constructing a causal hypothesis one should envisage as many *different* consequences of its truth as possible'' (p. 252). While rejecting a single point hypothesis will be satisfactory for some research projects, for others a set of findings is required to provide strong support for a theory. Predicting multiple outcomes provides a more stringent evaluation of the theory and is good scientific practice for these projects.

Traditional power analysis tools are designed for single hypothesis studies, but as noted, many scientific studies evaluate multiple hypotheses. The R package, **powermh**, estimates the sample size needed for studies with multiple hypotheses. Users go through four steps:

1. Define criteria for `success,'
2. Define data and the minimum effect sizes to detect,
3. Pose power analysis questions, and
4. Produce output, perhaps in the form of a plot.

The initial function (`pwAnova`) allows this for Oneway ANOVA designs. The user specifies the criteria for success (e.g., having all p-values being significant, having one p-value being significant but all in the predict direction), how the data are distributed and the minimum effect sizes to detect (e.g., log-normal and means .5sd higher than control), what parameter is to be output (e.g., the suggested sample size), and some limited control over any plots that can be created (the user may wish to create plots from the numeric matrix that is output). The package is at https://github.com/dbrookswr/powermh. Details of the pwAnova function are there, and functions for other designs are in development.

In summary, traditional power analysis reinforces the notion that the primary goal of a study is to achieve a single *p*-value less than some value thereby warranting publication. This focus on the significance of a single test statistic simplifies what constitutes `success' for many studies. Some designs can have only a single hypothesis to investigate. Where appropriate simple designs are preferred. Cohen's students describe this in the extreme: ``some of my students have spread the rumor that my idea of the perfect study is one with 10,000 cases and no variables. They go too far'' [-@Cohen1990, p.~1305]. This tongue-in-cheek ideal is in opposition to Fisher's quotation at the start of this paper and what can occur when considering many potential causes that scientists investigate. For many scientific projects it is appropriate to investigate a multitude of effects that may be caused by any manipulation. While power packages and tables exist for estimating power for simple one-hypothesis studies (so not quite as extreme as Cohen's students jokingly suggest), the **powermh** packages allows scientists to calculate power for studies with multiple hypotheses.


##Example

Further examples and details of the options are available: https://github.com/dbrookswr/powermh/blob/master/pwAnova.pdf. This is Example 2 from that document.

The researcher is asking if there are any significant pairwise differences for the means from a five-group Oneway ANOVA. As there are five conditions, there are 10 pairwise comparisons. This cannot be represented in the contrast matrix with 4 contrasts, so the `extrasuccess` slot is used in conjunction with `pairwise.t.test` function (from the **stats** package). Holm's method is used to adjust for the number of $p$ values calculated (this is the default for `pairwise.t.test`), so that the Type 1 error rate should be near the nominal rate. Holm's methods is preferred over competitors like Bonferroni's method because it maintains family-wise Type 1 error rate while having more power. The results are shown in the plots on the following page. For the left panels the null hypothesis is true. In the right panels $\mu_1 = 0; \mu_2 = \mu_3 = .1; \mu_4 = \mu_5 = .3$ in standard deviation units. All distributions are normal. The left panels show that Type 1 error is maintained, and right panels show the power for different sample sizes. The bottom row uses no adjustment for the number of *p*-values and the proportion of `successes' increases substantially.

```
reps <- 1000
dd3 <- list(function(n=nv,ef=efv,...) rnorm(n,ef),
            function(n=nv,ef=efv,...) rnorm(n,ef),
            function(n=nv,ef=efv,...) rnorm(n,ef),
            function(n=nv,ef=efv,...) rnorm(n,ef),
            function(n=nv,ef=efv,...) rnorm(n,ef))
par(mfrow=c(2,2))
eg3a <- pwAnova(dd3,dfv=4,
        replics=reps,varyef=rep(0,5),varyn=seq(100,500,20),r2p=.05,
        extrasuccess = list(
          function() any(pairwise.t.test(dd[,2],dd[,1],
              p.adjust.method="holm")$p.value < .05,na.rm=TRUE)))
eg3b <- pwAnova(dd3,dfv=4,
        replics=reps,varyef=c(0,.1,.1,.3,.3),varyn=seq(100,500,20),
        r2p=.05,extrasuccess = list(
          function() any(pairwise.t.test(dd[,2],dd[,1],
              p.adjust.method="holm")$p.value < .05,na.rm=TRUE)))
eg3c <- pwAnova(dd3,dfv=4,
        replics=reps,varyef=rep(0,5),varyn=seq(100,500,20),
        extrasuccess = list(
          function() any(pairwise.t.test(dd[,2],dd[,1],
              p.adjust.method="none")$p.value < .05,na.rm=TRUE)))
eg3d <- pwAnova(dd3,dfv=4,
        replics=reps,varyef=c(0,.1,.1,.3,.3),varyn=seq(100,500,20),
        extrasuccess = list( 
          function() any(pairwise.t.test(dd[,2],dd[,1],
              p.adjust.method="none")$p.value < .05,na.rm=TRUE)))
```
\pagebreak
The following plot is produced.

![Power Plots for Example 2.](eg2pwAnova.pdf)

# References

