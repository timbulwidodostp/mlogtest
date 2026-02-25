{smcl}
{* 2014-08-06 scott long & jeremy freese}{...}
{title:Title}

{p2colset 5 17 18 2}{...}
{p2col:{cmd:mlogtest} {hline 2}}Tests for {cmd:mlogit}{p_end}
{p2colreset}{...}


{title:General syntax}

{p 4 18 2}
{cmd:mlogtest} [{it:varlist}] [, {opt all} ... ]

{p 2 4 2}
{it:Options for tests of variables}
{p_end}

{p 6 18 2}
[ {opt w:ald} {opt lr} {opt set(setname1 = varlist1  \ setname2=varlist2 \ [...])} ]
{p_end}

{p 2 4 2}
{it:Options for tests of IIA}
{p_end}

{p 6 18 2}
[ {opt iia} {opt haus:man} {opt sm:hsiao} {opt suest} ]
{p_end}

{p 2 4 2}
{it:Options for tests for combining outcome}
{p_end}

{p 6 18 2}
[ {opt comb:ine} {opt lrcomb:ine} ]
{p_end}


{title:Overview}

{pstd}
{cmd:mlogtest} computes three types of tests for the multinomial logit model.
First, for each regressor it can perform an LR or
Wald test of the hypothesis that all coefficients for each regressor
are zero. Sets of variables can be tested simultaneously.
Second, it computes the Hausman-McFadden test of the
independence of irrelevance alternatives (IIA), the Small-Hsiao test, or
the SUEST version of the Hausman-McFadden test.
Third, it computes Wald and LR ratio tests of whether pair of outcome
categories can be combined.
Multiple tests
are available.

{pstd}
Output is much clearer if you assign value labels to the dependent variable.
See {help label}
{p_end}


{title:Options}

{p 4 4 2}
{it:Tests of variables}
{p_end}
{p2colset 8 20 20 0}
{synopt:{opt all}}All tests should be computed.

{synopt:{opt lr}}Compute LR tests for each regressor.

{synopt:{opt w:ald}}Compute Wald tests for each regressors.
{p_end}

{synopt:{opt set(setname1 = varset1 \ setname2 = varset2 ...)}}Specifies a set of
regressors is to be jointly tested with the
{cmd:lr} or {cmd:wald} options. The slash {cmd:\} divides
different sets.  This option can be used to test if coefficients for age and age-squared
are simultaneously equal to zero, or that coefficients for all indicators of a factor
variables are simultaneously equal to zero.
{p_end}

{p 4 4 2}
{it:Tests of IIA}
{p_end}

{synopt:{opt haus:man}}Hausman-McFadden tests using the {help hausman} command.
{p_end}

{synopt:{opt suest}}Hausman tests using the {help suest} command.
{p_end}

{synopt:{opt sm:hsiao}}Small-Hsiao tests. To reproduce results,
you should {cmd:set seed} before running this test.
{p_end}

{synopt:{opt iia}}Run all IIA tests.
{p_end}

{synopt:{opt detail}}Reports detailed output of the steps used to compute
the tests.
{p_end}

{p 4 4 2}
{it:Tests for combining outcome}
{p_end}

{synopt:{opt comb:ine}}Wald tests of whether two outcomes can be combined.
{p_end}

{synopt:{opt lrcomb:ine}}LR tests of whether two outcomes can be combined.
{p_end}


{title:Returned matrices}

{synopt:{opt r(combine)}}results of Wald tests to combine categories.
{p_end}

{synopt:{opt r(lrcomb)}}results of LR tests to combine categories.
{p_end}

{synopt:{opt r(hausman)}}results of Hausman tests of IIA assumption.
{p_end}

{synopt:{opt r(smhsiao)}}results of Small-Hsiao tests of IIA.
{p_end}

{synopt:{opt r(wald)}}results of Wald test that all coefficients for a regressor
are zero.
{p_end}

{synopt:{opt r(lrtest)}}results of LR test that all coefficients for a regressor
are zero.
{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
{ul:Example 1: Compute lr test for independent variables}{p_end}

{phang2}{cmd:. spex mlogit}{p_end}
{phang2}{cmd:. mlogtest, lr set(educ_set = 2.educ 3.educ)}{p_end}

{pstd}
{ul:Example 2: Test iia following mlogit}{p_end}

{phang2}{cmd:. spex mlogit}{p_end}
{phang2}{cmd:. set seed 85193017}{p_end}
{phang2}{cmd:. mlogtest, iia}{p_end}

{pstd}
{ul:Example 3: Test if outcome categories can be combined}{p_end}

{phang2}{cmd:. spex mlogit}{p_end}
{phang2}{cmd:. mlogtest, combine}{p_end}


{title:Acknowledgment}

{p 4 2 2}
Code used for the the Small-Hsiao test is based on a program by Nick Winter.

INCLUDE help spost13_footer
