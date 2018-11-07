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
     orcid: 0000-0003-0258-2158 
     affiliation: 2
affiliations:
   - name: Alder Graduate School of Education
     index: 1
   - name: Loyola University Maryland
     index: 2 
date: 5 July 2018
bibliography: paper.bib
---


# Summary


Power analysis is recommended by scientific societies and journals [@WilkinsonEA1999; @Wright2003], and required by some funding bodies (e.g., ies.ed.gov/funding/resources.asp). Tables [@Cohen1992] and computer packages [@FaulEA2007] can assist researchers in conducting traditional power analyses. The algorithms used for traditional power analysis require the researcher to choose a single test statistic.  

The reality of many research projects is more complex. Most processes in nature cause multiple effects, and as such a useful exploration of this process should make multiple predictions. Fisher summed this up nicely: "Make your theories elaborate" [@Cochran1965, p. 252]. Cochran explains this statement: "What Sir Ronald meant ... was that when constructing a causal hypothesis one should envisage as many *different* consequences of its truth as possible" (p. 252, *emphasis* in original). While rejecting a single point hypothesis will be satisfactory for some research projects, for others a set of findings is required to provide strong support for a theory. Predicting multiple outcomes provides a more stringent evaluation of the theory and is good scientific practice for these projects.

Traditional power analysis tools are designed for single hypothesis studies, but as noted, many scientific studies evaluate multiple hypotheses. The R package **powermh** estimates the sample size needed for studies with multiple hypotheses. Users go through four steps:

1. Define the criteria for success,
2. Define the data and the minimum effect sizes the study is designed to have a high likelihood to detect,
3. Pose power analysis questions, and
4. Produce output, perhaps in the form of a plot.

The function (`pwAnova`) allows this for ANOVA designs. The user specifies the criteria for success (e.g., having all p-values being significant, having one p-value being significant but all others in the predicted direction), how the data are distributed and the minimum effect sizes to detect (e.g., log-normal and means .5sd higher than control), what parameter is to be output (e.g., the suggested sample size), and some limited control over any plots that can be created (the user may wish to create plots from the numeric matrix that is output). The package is at https://github.com/dbrookswr/powermh and will be placed on CRAN. 

install_github("dbrookswr/powermh")

Details of the pwAnova function are there along with several examples. 

In summary, traditional power analysis reinforces the notion that the primary goal of a study is to achieve a single *p*-value less than some value thereby warranting publication. This focus on the significance of a single test statistic simplifies what constitutes *success* for many studies. Some designs can have only a single hypothesis to investigate. Where appropriate simple designs are preferred. Cohen's students describe this in the extreme: ``some of my students have spread the rumor that my idea of the perfect study is one with 10,000 cases and no variables. They go too far'' [-@Cohen1990, p.~1305]. This tongue-in-cheek ideal is in opposition to Fisher's quotation at the start of this paper and what can occur when considering many potential causes that scientists investigate. For many scientific projects it is appropriate to investigate a multitude of effects that may be caused by any manipulation. While power packages and tables exist for estimating power for simple one-hypothesis studies (so not quite as extreme as Cohen's students jokingly suggest), the **powermh** packages allows scientists to calculate power for studies with multiple hypotheses.


##Simple Example

Further (and more complex) examples and details of the options are available: https://github.com/dbrookswr/powermh/blob/master/pwAnova.pdf. This is Example 5 from that document.

The researcher has a sample of 500 and five independent groups. The researcher is interested if all four contrasts are significant at the default five percent level (pcon in the function below can be a vector with length equal to the number of hypotheses, or a single value). The default R contrasts are used (the first group is compared with the others as if the first is a control group). The varyn option tells the function how to vary the sample size. The values reported in the plot are for one condition.

First, the data are constructed. Here the data are normally distributed and each experimental group is half a standard deviation higher than the control group. The researcher is interested in how many people to have in each group. A more powerful design results from having a larger control group and this is shown in Example 5 from the pdf document listed above. A horizontal line is added to the plot at .8, the convention many aim for regarding power. This shows that you should have about 170 people in each condition.

```
options(warn=-1)
set.seed(223)
gr <- sample(1:5,500,replace=TRUE)
yvar <- rnorm(500) + .5 * as.numeric(gr != 1)
sampdata <- cbind(gr,yvar)
options(warn=-1)
eg <- pwAnova(sampdata,
   replics=1000,varyn=seq(675,1075,25),pcon=.05,dfv=3)
abline(h=.8,lty=2)
```

The following plot is produced.

![Power Plot for Example.](eg5.pdf)

If the research only requires one contrast to be significant a much smaller sample is suggested (about 28 per condition for 80 percent power). dfv control the smoothness of the curve (by controlling the degrees of freedom of the spline) in the plot. 

```
egB <- pwAnova(sampdata,
   replics=1000,varyn=seq(100,175,5),pcon=.05,numcon=1,dfv=3)
abline(h=.8,lty=2)
```

![Power Plot for Example, part B.](eg5b.pdf)


# References

