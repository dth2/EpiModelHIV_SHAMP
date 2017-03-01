
# SHAMP  -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param last.neg.test.B.int Time range in days for last negative test for
#'        black men.
#' @param mean.test.B.int Mean intertest interval in days for black MSM who test.
#' @param last.neg.test.W.int Time range in days for last negative test for
#'        white men.
#' @param mean.test.W.int Mean intertest interval in days for white MSM who test.
#' @param testing.pattern Method for HIV testing, with options \code{"memoryless"}
#'        for constant hazard without regard to time since previous test, or
#'        \code{"interval"} deterministic fixed intervals.
#' @param test.window.int Length of the HIV test window period in days.
#' @param tt.traj.B.prob Proportion of black MSM who enter one of four
#'        testing/treatment trajectories: never test or treat, test and never
#'        initiate treatment, test and treated with partial viral suppression,
#'        and test and treated with full suppression.
#' @param tt.traj.W.prob Proportion of white MSM who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tx.init.B.prob Probability per time step that a black MSM who has
#'        tested positive will initiate treatment.
#' @param tx.init.W.prob Probability per time step that a white MSM who has
#'        tested positive will initiate treatment.
#' @param tx.halt.B.prob Probability per time step that a black MSM who is
#'        currently on treatment will halt treatment.
#' @param tx.halt.W.prob Probability per time step that a white MSM who is
#'        currently on treatment will halt treatment.
#' @param tx.reinit.B.prob Probability per time step that a black MSM who is
#'        not currently on treatment but who has been in the past will
#'        re-initiate treatment.
#' @param tx.reinit.W.prob Probability per time step that a white MSM who is
#'        not currently on treatment but who has been in the past will
#'        re-initiate treatment.
#' @param max.time.off.tx.full.int Number of days off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.part.int Number of days on treatment for a
#'        partial suppressor beofre onset of AIDS.
#' @param max.time.off.tx.part.int Nnumber of days off treatment for a
#'        partial suppressor before onset of AIDS, including time before
#'        diagnosis.
#' @param vl.acute.rise.int Number of days to peak viremia during acute
#'        infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute
#'        infection.
#' @param vl.acute.fall.int Number of days from peak viremia to set-point
#'        viral load during the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of days to AIDS for a treatment-naive
#'        patient.
#' @param vl.aids.int Duration of AIDS stage infection in days.
#' @param vl.fatal Viral load in AIDS at which death occurs.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param full.supp.down.slope For full suppressors, number of log10 units that
#'        viral load falls per time step from treatment initiation or re-initiation
#'        until the level in \code{vl.full.supp}.
#' @param full.supp.up.slope For full suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected
#'        value.
#' @param part.supp.down.slope For partial suppressors, number of log10 units
#'        that viral load falls per time step from treatment initiation or
#'        re-initiation until the level in \code{vl.part.supp}.
#' @param part.supp.up.slope For partial suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected value.
#' @param b.B.rate Rate at which black MSM enter the population.
#' @param b.W.rate Rate at which white MSM enter the population.
#' @param birth.age Age (in years) of new arrivals.
#' @param b.method Method for calculating the number of expected births at each
#'        time step, with \code{"fixed"} based on the number of persons at the
#'        initial time step and \code{"varying"} based on the current time step.
#' @param URAI.prob Probability of transmission for a man having unprotected
#'        receptive anal intercourse with an infected man at set point viral
#'        load.
#' @param UIAI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive anal intercourse with an infected man at set
#'        point viral load.
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#' @param condom.rr Relative risk of infection from anal sex when a condom is
#'        used.
#' @param disc.outset.main.B.prob Probability that an HIV-infected black MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.outset.main.W.prob Probability that an HIV-infected white MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.at.diag.main.B.prob Probability that a black MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.at.diag.main.W.prob Probability that a white MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.main.B.prob Probability that an HIV-infected black MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.post.diag.main.W.prob Probability that an HIV-infected white MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.outset.pers.B.prob Probability that an HIV-infected black MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.outset.pers.W.prob Probability that an HIV-infected white MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.at.diag.pers.B.prob Probability that a black MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.at.diag.pers.W.prob Probability that a white MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.pers.B.prob Probability that an HIV-infected black MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.post.diag.pers.W.prob Probability that an HIV-infected white MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.inst.B.prob Probability that an HIV-infected black MSM will
#'        disclose his status to a one-off partner.
#' @param disc.inst.W.prob Probability that an HIV-infected white MSM will
#'        disclose his status to a one-off partner.
#' @param circ.B.prob Probablity that a black new arrival in the population
#'        will be circumcised.
#' @param circ.W.prob Probablity that a white new arrival in the population
#'        will be circumcised.
#' @param ccr5.B.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among black MSM.
#' @param ccr5.W.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among white MSM.
#' @param ccr5.heteroz.rr Relative risk of infection for men who are heterozygous
#'        in the CCR5 mutation.
#' @param num.inst.ai.classes Number of quantiles into which men should be
#'        divided in determining their levels of one-off anal intercourse.
#' @param base.ai.main.BB.rate Expected coital frequency in black-black main
#'        partnerships (acts per day).
#' @param base.ai.main.BW.rate Expected coital frequency in black-white main
#'        partnerships (acts per day).
#' @param base.ai.main.WW.rate Expected coital frequency in white-white main
#'        partnerships (acts per day).
#' @param base.ai.pers.BB.rate Expected coital frequency in black-black casual
#'        partnerships (acts per day).
#' @param base.ai.pers.BW.rate Expected coital frequency in black-white casual
#'        partnerships (acts per day).
#' @param base.ai.pers.WW.rate Expected coital frequency in white-white casual
#'        partnerships (acts per day).
#' @param ai.scale General relative scaler for all act rates for model
#'        calibration.
#' @param cond.main.BB.prob Probability of condom use in a black-black main
#'        partnership.
#' @param cond.main.BW.prob Probability of condom use in a black-white main
#'        partnership.
#' @param cond.main.WW.prob Probability of condom use in a white-white main
#'        partnership.
#' @param cond.pers.always.prob Fraction of men in casual partnerships who always
#'        use condoms in those partnerships.
#' @param cond.pers.BB.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-black casual partnerships.
#' @param cond.pers.BW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-white casual partnerships.
#' @param cond.pers.WW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a white-white casual partnerships.
#' @param cond.inst.always.prob Fraction of men in instant partnerships who always
#'        use condoms in those partnerships.
#' @param cond.inst.BB.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-black one-off partnerships.
#' @param cond.inst.BW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-white one-off partnerships.
#' @param cond.inst.WW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a white-white one-off partnerships.
#' @param cond.always.prob.corr Correlation coefficient for probability of always
#'        using condoms in both casual and one-off
#' @param cond.rr.BB Condom probability scaler for black-black partnerships for
#'        model calibration purposes.
#' @param cond.rr.BW Condom probability scaler for black-white partnerships for
#'        model calibration purposes.
#' @param cond.rr.WW Condom probability scaler for white-white partnerships for
#'        model calibration purposes.
#' @param cond.diag.main.beta Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.main.beta Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV-infected man has disclosed.
#' @param cond.diag.pers.beta Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.pers.beta Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV-infected man has disclosed
#'        his status.
#' @param cond.diag.inst.beta Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.inst.beta Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV-infected man has disclosed
#'        his status.
#' @param vv.iev.BB.prob Probability that in a black-black partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param vv.iev.BW.prob Probability that in a black-white partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param vv.iev.WW.prob Probability that in a white-white partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.elig.model Modeling approach for determining who is eligible for
#'        PrEP. Current options are limited to: \code{"all"} for all persons who
#'        have never been on PrEP and are disease-susceptible.
#' @param prep.class.prob The probability of adherence class in non-adherent,
#'        low adherence, medium adherence, or high adherence groups (from Liu).
#' @param prep.class.hr The hazard ratio for infection per act associated with each
#'        level of adherence (from Grant).
#' @param prep.coverage The proportion of the eligible population who are start
#'        PrEP once they become eligible.
#' @param prep.cov.method The method for calculating PrEP coverage, with options
#'        of \code{"curr"} to base the numerator on the number of people currently
#'        on PrEP and \code{"ever"} to base it on the number of people ever on
#'        PrEP.
#' @param prep.cov.rate The rate at which persons initiate PrEP conditional on
#'        their eligibility, with 1 equal to instant start.
#' @param prep.tst.int Testing interval for those who are actively on PrEP. This
#'        overrides the mean testing interval parameters.
#' @param prep.risk.int Time window for assessment of risk eligibility for PrEP
#'        in days.
#' @param prep.risk.reassess If \code{TRUE}, reassess eligibility for PrEP at
#'        each testing visit.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param_shamp}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords shamp
#'
#' @export
param_shamp <- function(race.method = 1,
                      time.unit = 7,
                      last.neg.test.B.f.int = 301,
                      last.neg.test.BI.f.int = 301,
                      last.neg.test.H.f.int = 301,
                      last.neg.test.HI.f.int = 301,
                      last.neg.test.W.f.int = 30,
                      last.neg.test.B.msf.int = 301,
                      last.neg.test.BI.msf.int = 301,
                      last.neg.test.H.msf.int = 301,
                      last.neg.test.HI.msf.int = 301,
                      last.neg.test.W.msf.int = 30,
                      last.neg.test.B.msm.int = 301,
                      last.neg.test.BI.msm.int = 301,
                      last.neg.test.H.msm.int = 301,
                      last.neg.test.HI.msm.int = 301,
                      last.neg.test.W.msm.int = 301,
                      last.neg.test.B.msmf.int = 301,
                      last.neg.test.BI.msmf.int = 301,
                      last.neg.test.H.msmf.int = 301,
                      last.neg.test.HI.msmf.int = 301,
                      last.neg.test.W.msmf.int = 301,
                      mean.test.B.f.int = 301,
                      mean.test.BI.f.int = 301,
                      mean.test.H.f.int = 301,
                      mean.test.HI.f.int = 301,
                      mean.test.W.f.int = 30,
                      mean.test.B.msf.int = 301,
                      mean.test.BI.msf.int = 301,
                      mean.test.H.msf.int = 301,
                      mean.test.HI.msf.int = 301,
                      mean.test.W.msf.int = 30,
                      mean.test.B.msm.int = 301,
                      mean.test.BI.msm.int = 301,
                      mean.test.H.msm.int = 301,
                      mean.test.HI.msm.int = 301,
                      mean.test.W.msm.int = 301,
                      mean.test.B.msmf.int = 301,
                      mean.test.BI.msmf.int = 301,
                      mean.test.H.msmf.int = 301,
                      mean.test.HI.msmf.int = 301,
                      mean.test.W.msmf.int = 301,
                      testing.pattern = "memoryless",
                      test.window.int = 21,

                      tt.traj.B.f.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.BI.f.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.H.f.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.HI.f.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.W.f.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.B.msf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.BI.msf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.H.msf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.HI.msf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.W.msf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.B.msm.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.BI.msm.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.H.msm.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.HI.msm.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.W.msm.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.B.msmf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.BI.msmf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.H.msmf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.HI.msmf.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.W.msmf.prob = c(0.077, 0.000, 0.356, 0.567),
                      
                      tx.init.B.f.prob = 0.092,
                      tx.init.BI.f.prob = 0.092,
                      tx.init.H.f.prob = 0.092,
                      tx.init.HI.f.prob = 0.092,
                      tx.init.W.f.prob = 0.092,
                      tx.init.B.msf.prob = 0.092,
                      tx.init.BI.msf.prob = 0.092,
                      tx.init.H.msf.prob = 0.092,
                      tx.init.HI.msf.prob = 0.092,
                      tx.init.W.msf.prob = 0.092,
                      tx.init.B.msm.prob = 0.092,
                      tx.init.BI.msm.prob = 0.092,
                      tx.init.H.msm.prob = 0.092,
                      tx.init.HI.msm.prob = 0.092,
                      tx.init.W.msm.prob = 0.092,
                      tx.init.B.msmf.prob = 0.092,
                      tx.init.BI.msmf.prob = 0.092,
                      tx.init.H.msmf.prob = 0.092,
                      tx.init.HI.msmf.prob = 0.092,
                      tx.init.W.msmf.prob = 0.092,

                      tx.halt.B.f.prob = 0.0102,
                      tx.halt.BI.f.prob = 0.0102,
                      tx.halt.H.f.prob = 0.0102,
                      tx.halt.HI.f.prob = 0.0102,
                      tx.halt.W.f.prob = 0.0102,
                      tx.halt.B.msf.prob = 0.0102,
                      tx.halt.BI.msf.prob = 0.0102,
                      tx.halt.H.msf.prob = 0.0102,
                      tx.halt.HI.msf.prob = 0.0102,
                      tx.halt.W.msf.prob = 0.0102,
                      tx.halt.B.msm.prob = 0.0102,
                      tx.halt.BI.msm.prob = 0.0102,
                      tx.halt.H.msm.prob = 0.0102,
                      tx.halt.HI.msm.prob = 0.0102,
                      tx.halt.W.msm.prob = 0.0102,
                      tx.halt.B.msmf.prob = 0.0102,
                      tx.halt.BI.msmf.prob = 0.0102,
                      tx.halt.H.msmf.prob = 0.0102,
                      tx.halt.HI.msmf.prob = 0.0102,
                      tx.halt.W.msmf.prob = 0.0102,

                      tx.reinit.B.f.prob = 0.00066,
                      tx.reinit.BI.f.prob = 0.00066,
                      tx.reinit.H.f.prob = 0.00066,
                      tx.reinit.HI.f.prob = 0.00066,
                      tx.reinit.W.f.prob = 0.00066,
                      tx.reinit.B.msf.prob = 0.00066,
                      tx.reinit.BI.msf.prob = 0.00066,
                      tx.reinit.H.msf.prob = 0.00066,
                      tx.reinit.HI.msf.prob = 0.00066,
                      tx.reinit.W.msf.prob = 0.00066,
                      tx.reinit.B.msm.prob = 0.00066,
                      tx.reinit.BI.msm.prob = 0.00066,
                      tx.reinit.H.msm.prob = 0.00066,
                      tx.reinit.HI.msm.prob = 0.00066,
                      tx.reinit.W.msm.prob = 0.00066,
                      tx.reinit.B.msmf.prob = 0.00066,
                      tx.reinit.BI.msmf.prob = 0.00066,
                      tx.reinit.H.msmf.prob = 0.00066,
                      tx.reinit.HI.msmf.prob = 0.00066,
                      tx.reinit.W.msmf.prob = 0.00066,

                      max.time.off.tx.full.int = 520 * 7,
                      max.time.on.tx.part.int = 52 * 15 * 7,
                      max.time.off.tx.part.int = 520 * 7,
                      vl.acute.rise.int = 45,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 45,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520 * 7,
                      vl.aids.int = 52 * 2 * 7,
                      vl.fatal = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      full.supp.down.slope = 0.25,
                      full.supp.up.slope = 0.25,
                      part.supp.down.slope = 0.25,
                      part.supp.up.slope = 0.25,

                      b.B.f.rate = 1e-3 / 7,
                      b.BI.f.rate = 1e-3 / 7,
                      b.H.f.rate = 1e-3 / 7,
                      b.HI.f.rate = 1e-3 / 7,
                      b.W.f.rate = 1e-3 / 7,
                      b.B.m.rate = 1e-3 / 7,
                      b.BI.m.rate = 1e-3 / 7,
                      b.H.m.rate = 1e-3 / 7,
                      b.HI.m.rate = 1e-3 / 7,
                      b.W.m.rate = 1e-3 / 7,

                      birth.age = 18,
                      b.method = "fixed",
                      msm.frac=.03,
                      msmf.frac=.1,

                      URAI.prob = 0.0082 * 1.09,
                      UIAI.prob = 0.0031 * 1.09,
                      
                      URVI.prob = 0.0082 * 1.09,
                      UIVI.prob = 0.0031 * 1.09,
                      
                      acute.rr = 6,
                      circ.rr = 0.4,
                      condom.rr = 0.295,

                      disc.outset.main.B.f.prob = 0.685,
                      disc.outset.main.BI.f.prob = 0.685,
                      disc.outset.main.H.f.prob = 0.685,
                      disc.outset.main.HI.f.prob = 0.685,
                      disc.outset.main.W.f.prob = 0.685,
                      disc.outset.main.B.msf.prob = 0.685,
                      disc.outset.main.BI.msf.prob = 0.685,
                      disc.outset.main.H.msf.prob = 0.685,
                      disc.outset.main.HI.msf.prob = 0.685,
                      disc.outset.main.W.msf.prob = 0.685,
                      disc.outset.main.B.msm.prob = 0.685,
                      disc.outset.main.BI.msm.prob = 0.685,
                      disc.outset.main.H.msm.prob = 0.685,
                      disc.outset.main.HI.msm.prob = 0.685,
                      disc.outset.main.W.msm.prob = 0.685,
                      disc.outset.main.B.msmf.prob = 0.685,
                      disc.outset.main.BI.msmf.prob = 0.685,
                      disc.outset.main.H.msmf.prob = 0.685,
                      disc.outset.main.HI.msmf.prob = 0.685,
                      disc.outset.main.W.msmf.prob = 0.685,
                      
                      disc.at.diag.main.B.f.prob = 1,
                      disc.at.diag.main.BI.f.prob = 1,
                      disc.at.diag.main.H.f.prob = 1,
                      disc.at.diag.main.HI.f.prob = 1,
                      disc.at.diag.main.W.f.prob = 1,
                      disc.at.diag.main.B.msf.prob = 1,
                      disc.at.diag.main.BI.msf.prob = 1,
                      disc.at.diag.main.H.msf.prob = 1,
                      disc.at.diag.main.HI.msf.prob = 1,
                      disc.at.diag.main.W.msf.prob = 1,
                      disc.at.diag.main.B.msm.prob = 1,
                      disc.at.diag.main.BI.msm.prob = 1,
                      disc.at.diag.main.H.msm.prob = 1,
                      disc.at.diag.main.HI.msm.prob = 1,
                      disc.at.diag.main.W.msm.prob = 1,
                      disc.at.diag.main.B.msmf.prob = 1,
                      disc.at.diag.main.BI.msmf.prob = 1,
                      disc.at.diag.main.H.msmf.prob = 1,
                      disc.at.diag.main.HI.msmf.prob = 1,
                      disc.at.diag.main.W.msmf.prob = 1,
                      
                      disc.post.diag.main.B.f.prob = 0,
                      disc.post.diag.main.BI.f.prob = 0,
                      disc.post.diag.main.H.f.prob = 0,
                      disc.post.diag.main.HI.f.prob = 0,
                      disc.post.diag.main.W.f.prob = 0,
                      disc.post.diag.main.B.msf.prob = 0,
                      disc.post.diag.main.BI.msf.prob = 0,
                      disc.post.diag.main.H.msf.prob = 0,
                      disc.post.diag.main.HI.msf.prob = 0,
                      disc.post.diag.main.W.msf.prob = 0,
                      disc.post.diag.main.B.msm.prob = 0,
                      disc.post.diag.main.BI.msm.prob = 0,
                      disc.post.diag.main.H.msm.prob = 0,
                      disc.post.diag.main.HI.msm.prob = 0,
                      disc.post.diag.main.W.msm.prob = 0,
                      disc.post.diag.main.B.msmf.prob = 0,
                      disc.post.diag.main.BI.msmf.prob = 0,
                      disc.post.diag.main.H.msmf.prob = 0,
                      disc.post.diag.main.HI.msmf.prob = 0,
                      disc.post.diag.main.W.msmf.prob = 0,
                      
                      
                      disc.outset.pers.B.f.prob = 0.527,
                      disc.outset.pers.BI.f.prob = 0.527,
                      disc.outset.pers.H.f.prob = 0.527,
                      disc.outset.pers.HI.f.prob = 0.527,
                      disc.outset.pers.W.f.prob = 0.828,
                      disc.outset.pers.B.msf.prob = 0.527,
                      disc.outset.pers.BI.msf.prob = 0.527,
                      disc.outset.pers.H.msf.prob = 0.527,
                      disc.outset.pers.HI.msf.prob = 0.527,
                      disc.outset.pers.W.msf.prob = 0.828,
                      disc.outset.pers.B.msm.prob = 0.527,
                      disc.outset.pers.BI.msm.prob = 0.527,
                      disc.outset.pers.H.msm.prob = 0.527,
                      disc.outset.pers.HI.msm.prob = 0.527,
                      disc.outset.pers.W.msm.prob = 0.828,
                      disc.outset.pers.B.msmf.prob = 0.527,
                      disc.outset.pers.BI.msmf.prob = 0.527,
                      disc.outset.pers.H.msmf.prob = 0.527,
                      disc.outset.pers.HI.msmf.prob = 0.527,
                      disc.outset.pers.W.msmf.prob = 0.828,
                      
                      disc.at.diag.pers.B.f.prob = 1,
                      disc.at.diag.pers.BI.f.prob = 1,
                      disc.at.diag.pers.H.f.prob = 1,
                      disc.at.diag.pers.HI.f.prob = 1,
                      disc.at.diag.pers.W.f.prob = 1,
                      disc.at.diag.pers.B.msf.prob = 1,
                      disc.at.diag.pers.BI.msf.prob = 1,
                      disc.at.diag.pers.H.msf.prob = 1,
                      disc.at.diag.pers.HI.msf.prob = 1,
                      disc.at.diag.pers.W.msf.prob = 1,
                      disc.at.diag.pers.B.msm.prob = 1,
                      disc.at.diag.pers.BI.msm.prob = 1,
                      disc.at.diag.pers.H.msm.prob = 1,
                      disc.at.diag.pers.HI.msm.prob = 1,
                      disc.at.diag.pers.W.msm.prob = 1,
                      disc.at.diag.pers.B.msmf.prob = 1,
                      disc.at.diag.pers.BI.msmf.prob = 1,
                      disc.at.diag.pers.H.msmf.prob = 1,
                      disc.at.diag.pers.HI.msmf.prob = 1,
                      disc.at.diag.pers.W.msmf.prob = 1,
                      
                      disc.post.diag.pers.B.f.prob = 0,
                      disc.post.diag.pers.BI.f.prob = 0,
                      disc.post.diag.pers.H.f.prob = 0,
                      disc.post.diag.pers.HI.f.prob = 0,
                      disc.post.diag.pers.W.f.prob = 0,
                      disc.post.diag.pers.B.mf.prob = 0,
                      disc.post.diag.pers.BI.mf.prob = 0,
                      disc.post.diag.pers.H.mf.prob = 0,
                      disc.post.diag.pers.HI.mf.prob = 0,
                      disc.post.diag.pers.W.mf.prob = 0, 
                      disc.post.diag.pers.B.msm.prob = 0,
                      disc.post.diag.pers.BI.msm.prob = 0,
                      disc.post.diag.pers.H.msm.prob = 0,
                      disc.post.diag.pers.HI.msm.prob = 0,
                      disc.post.diag.pers.W.msm.prob = 0, 
                      disc.post.diag.pers.B.msmf.prob = 0,
                      disc.post.diag.pers.BI.msmf.prob = 0,
                      disc.post.diag.pers.H.msmf.prob = 0,
                      disc.post.diag.pers.HI.msmf.prob = 0,
                      disc.post.diag.pers.W.msmf.prob = 0, 
                      
                      disc.inst.B.f.prob = 0,
                      disc.inst.BI.f.prob = 0,
                      disc.inst.H.f.prob = 0,
                      disc.inst.HI.f.prob = 0,
                      disc.inst.W.f.prob = 0,
                      disc.inst.B.msf.prob = 0,
                      disc.inst.BI.msf.prob = 0,
                      disc.inst.H.msf.prob = 0,
                      disc.inst.HI.msf.prob = 0,
                      disc.inst.W.msf.prob = 0,
                      disc.inst.B.msm.prob = 0,
                      disc.inst.BI.msm.prob = 0,
                      disc.inst.H.msm.prob = 0,
                      disc.inst.HI.msm.prob = 0,
                      disc.inst.W.msm.prob = 0,
                      disc.inst.B.msmf.prob = 0,
                      disc.inst.BI.msmf.prob = 0,
                      disc.inst.H.msmf.prob = 0,
                      disc.inst.HI.msmf.prob = 0,
                      disc.inst.W.msmf.prob = 0,

                      circ.B.prob = 0.874,
                      circ.BI.prob = 0.874,
                      circ.H.prob = 0.918,
                      circ.HI.prob = 0.918,
                      circ.W.prob = 0.918,
                      
                      ccr5.B.f.prob = c(0, 0.034),
                      ccr5.BI.f.prob = c(0, 0.034),
                      ccr5.H.f.prob = c(0, 0.034),
                      ccr5.HI.f.prob = c(0, 0.034),
                      ccr5.W.f.prob = c(0.021, 0.176),
                      ccr5.B.m.prob = c(0, 0.034),
                      ccr5.BI.m.prob = c(0, 0.034),
                      ccr5.H.m.prob = c(0, 0.034),
                      ccr5.HI.m.prob = c(0, 0.034),
                      ccr5.W.m.prob = c(0.021, 0.176),
                      ccr5.heteroz.rr = 0.3,

                      num.inst.ai.classes = 1,
                      base.ai.main.B.rate = 0.17,
                      base.ai.main.BI.rate = 0.26,
                      base.ai.main.H.rate = 0.23,
                      base.ai.main.HI.rate = 0.26,
                      base.ai.main.W.rate = 0.23,
                      base.ai.pers.B.rate = 0.11,
                      base.ai.pers.BI.rate = 0.16,
                      base.ai.pers.H.rate = 0.14,
                      base.ai.pers.HI.rate = 0.14,
                      base.ai.pers.W.rate = 0.14,
                      
                      base.vi.main.B.rate = 0.17,
                      base.vi.main.BI.rate = 0.26,
                      base.vi.main.H.rate = 0.23,
                      base.vi.main.HI.rate = 0.26,
                      base.vi.main.W.rate = 0.23,
                      base.vi.pers.B.rate = 0.11,
                      base.vi.pers.BI.rate = 0.16,
                      base.vi.pers.H.rate = 0.14,
                      base.vi.pers.HI.rate = 0.16,
                      base.vi.pers.W.rate = 0.14,
                      ai.scale = 1,
                      vi.scale = 1,

                      cond.main.B.prob.msm = 0.38,
                      cond.main.BI.prob.msm = 0.10,
                      cond.main.H.prob.msm = 0.15,
                      cond.main.HI.prob.msm = 0.15,
                      cond.main.W.prob.msm = 0.15,
                      cond.pers.always.prob.msm = 0.216,
                      cond.pers.B.prob.msm = 0.26,
                      cond.pers.BI.prob.msm = 0.26,
                      cond.pers.H.prob.msm = 0.26,
                      cond.pers.HI.prob.msm = 0.26,
                      cond.pers.W.prob.msm = 0.26,
                      cond.inst.always.prob.msm = 0.326,
                      cond.inst.B.prob.msm = 0.27,
                      cond.inst.BI.prob.msm = 0.27,
                      cond.inst.H.prob.msm = 0.27,
                      cond.inst.HI.prob.msm = 0.27,
                      cond.inst.W.prob.msm = 0.27,
                      cond.always.prob.corr.msm = 0.5,
                      
                      cond.main.B.prob.het = 0.38,
                      cond.main.BI.prob.het = 0.10,
                      cond.main.H.prob.het = 0.15,
                      cond.main.HI.prob.het = 0.15,
                      cond.main.W.prob.het = 0.15,
                      cond.pers.always.prob.het = 0.216,
                      cond.pers.B.prob.het = 0.26,
                      cond.pers.BI.prob.het = 0.26,
                      cond.pers.H.prob.het = 0.26,
                      cond.pers.HI.prob.het = 0.26,
                      cond.pers.W.prob.het = 0.26,
                      cond.inst.always.prob.het = 0.326,
                      cond.inst.B.prob.het = 0.27,
                      cond.inst.BI.prob.het = 0.27,
                      cond.inst.H.prob.het = 0.27,
                      cond.inst.HI.prob.het = 0.27,
                      cond.inst.W.prob.het = 0.27,
                      cond.always.prob.corr.het = 0.5,
                      
                      
                      cond.rr.BB = 1,
                      cond.rr.BW = 1,
                      cond.rr.WW = 1,
                      cond.diag.main.beta.msm = -0.67,
                      cond.discl.main.beta.msm = -0.85,
                      cond.diag.pers.beta.msm = -0.67,
                      cond.discl.pers.beta.msm = -0.85,  
                      cond.diag.inst.beta.msm = -0.67,
                      cond.discl.inst.beta.msm = -0.85,
                      
                      cond.diag.main.beta.het = -0.67,
                      cond.discl.main.beta.het = -0.85,
                      cond.diag.pers.beta.het = -0.67,
                      cond.discl.pers.beta.het = -0.85,  
                      cond.diag.inst.beta.het = -0.67,
                      cond.discl.inst.beta.het = -0.85,
                      
                      role.B.msm.prob = c(0.242, 0.321, 0.437),
                      role.BI.msm.prob = c(0.242, 0.321, 0.437),
                      role.H.msm.prob = c(0.242, 0.321, 0.437),
                      role.HI.msm.prob = c(0.242, 0.321, 0.437),
                      role.W.msm.prob = c(0.242, 0.321, 0.437),
                      role.B.msmf.prob = c(0.242, 0.321, 0.437),
                      role.BI.msmf.prob = c(0.242, 0.321, 0.437),
                      role.H.msmf.prob = c(0.242, 0.321, 0.437),
                      role.HI.msmf.prob = c(0.242, 0.321, 0.437),
                      role.W.msmf.prob = c(0.242, 0.321, 0.437),

                      vv.iev.B.prob = 0.42,
                      vv.iev.BI.prob = 0.56,
                      vv.iev.H.prob = 0.56,
                      vv.iev.HI.prob = 0.56,
                      vv.iev.W.prob = 0.49,

                      prep.start = Inf,
                      prep.elig.model = "base",
                      prep.class.prob = c(0.211, 0.07, 0.1, 0.619),
                      prep.class.hr = c(1, 0.69, 0.19, 0.05),
                      prep.coverage = 0,
                      prep.cov.method = "curr",
                      prep.cov.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 182,
                      prep.risk.reassess = TRUE,
                      
                      immig.depart.BI.f=0,
                      immig.depart.HI.f=0,
                      immig.depart.BI.m=0,
                      immig.depart.HI.m=0,
                      immig.return.BI.f=0,
                      immig.return.HI.f=0,
                      immig.return.BI.m=0,
                      immig.return.HI.m=0,
                      immig.aq.prob.BI.f=0,
                      immig.aq.prob.HI.f=0,
                      immig.aq.prob.BI.m=0,
                      immig.aq.prob.HI.m=0,  
                      
                      msm.aq.prob.B=0,
                      msm.aq.prob.BI=0,
                      msm.aq.prob.H=0,
                      msm.aq.prob.HI=0,
                      msm.aq.prob.W=0,

                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
          call. = FALSE)
  }


#  if (race.method == 1) {
#    p$last.neg.test.B.int = (last.neg.test.B.int + last.neg.test.W.int)/2
#    p$last.neg.test.W.int = (last.neg.test.B.int + last.neg.test.W.int)/2
#    p$mean.test.B.int = (mean.test.W.int + mean.test.B.int)/2
#    p$mean.test.W.int = (mean.test.W.int + mean.test.B.int)/2
#    p$tt.traj.B.prob = (tt.traj.B.prob + tt.traj.W.prob)/2
#    p$tt.traj.W.prob = (tt.traj.B.prob + tt.traj.W.prob)/2
#    p$tx.init.B.prob = (tx.init.B.prob + tx.init.W.prob)/2
#    p$tx.init.W.prob = (tx.init.B.prob + tx.init.W.prob)/2
#    p$tx.halt.B.prob = (tx.halt.B.prob + tx.halt.W.prob)/2
#    p$tx.halt.W.prob = (tx.halt.B.prob + tx.halt.W.prob)/2
#    p$tx.reinit.B.prob = (tx.reinit.B.prob + tx.reinit.W.prob)/2
#    p$tx.reinit.W.prob = (tx.reinit.B.prob + tx.reinit.W.prob)/2
#    p$disc.outset.main.B.prob = (disc.outset.main.B.prob + disc.outset.main.W.prob)/2
#    p$disc.outset.main.W.prob = (disc.outset.main.B.prob + disc.outset.main.W.prob)/2
#    p$disc.outset.pers.B.prob = (disc.outset.pers.B.prob + disc.outset.pers.W.prob)/2
#    p$disc.outset.pers.W.prob = (disc.outset.pers.B.prob + disc.outset.pers.W.prob)/2
#    p$disc.inst.B.prob = (disc.inst.B.prob + disc.inst.W.prob)/2
#    p$disc.inst.W.prob = (disc.inst.B.prob + disc.inst.W.prob)/2
#    p$circ.B.prob = (circ.B.prob + circ.W.prob)/2
#    p$circ.W.prob = (circ.B.prob + circ.W.prob)/2
#    p$ccr5.B.prob = (ccr5.B.prob + ccr5.W.prob)/2
#    p$ccr5.W.prob = (ccr5.B.prob + ccr5.W.prob)/2
#    p$base.ai.main.BB.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
#                                base.ai.main.WW.rate)/3
#    p$base.ai.main.BW.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
#                                base.ai.main.WW.rate)/3
#    p$base.ai.main.WW.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
#                                base.ai.main.WW.rate)/3
#    p$base.ai.pers.BB.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
#                                base.ai.pers.WW.rate)/3
#    p$base.ai.pers.BW.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
#                                base.ai.pers.WW.rate)/3
#    p$base.ai.pers.WW.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
#                                base.ai.pers.WW.rate)/3
#    p$cond.main.BB.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
#    p$cond.main.BW.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
#    p$cond.main.WW.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
#    p$cond.pers.BB.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
#    p$cond.pers.BW.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
#    p$cond.pers.WW.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
#    p$cond.inst.BB.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
#    p$cond.inst.BW.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
#    p$cond.inst.WW.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
#    p$vv.iev.BB.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
#    p$vv.iev.BW.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
#    p$vv.iev.WW.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
#  }

  p$time.unit <- time.unit

  intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
  p[intvars] <- lapply(p[intvars], FUN = function(x) round(x / p$time.unit))

  ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
  p[ratevars] <- lapply(p[ratevars], FUN = function(x) x * p$time.unit)

  
  p$role.B.msm.prob <- role.B.msm.prob
  p$role.BI.msm.prob <- role.BI.msm.prob
  p$role.H.msm.prob <- role.H.msm.prob
  p$role.HI.msm.prob <- role.HI.msm.prob
  p$role.W.msm.prob <- role.W.msm.prob
  
  p$role.B.msmf.prob <- role.B.msmf.prob
  p$role.BI.msmf.prob <- role.BI.msmf.prob
  p$role.H.msmf.prob <- role.H.msmf.prob
  p$role.HI.msmf.prob <- role.HI.msmf.prob
  p$role.W.msmf.prob <- role.W.msmf.prob

  p$inst.trans.matrix <- matrix(1, nrow = 1)
  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)


  p$riskh.start <- max(1, prep.start - prep.risk.int - 1)

  p$method <- method
  p$modes <- 1

  p$asmr.B.f <- asmr.B.f
  p$asmr.BI.f <- asmr.BI.f
  p$asmr.H.f <- asmr.H.f
  p$asmr.HI.f <- asmr.HI.f
  p$asmr.W.f <- asmr.W.f
  p$asmr.B.m <- asmr.B.m
  p$asmr.BI.m <- asmr.BI.m
  p$asmr.H.m <- asmr.H.m
  p$asmr.HI.m <- asmr.HI.m
  p$asmr.W.m <- asmr.W.m


  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param prev.B Initial disease prevalence among black MSM.
#' @param prev.W Initial disease prevalence among white MSM.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init_shamp}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @keywords SHAMP
#'
#' @export
init_shamp <- function(prev.B.f = 0.005,
                     prev.BI.f =0.005,
                     prev.H.f =0.005,
                     prev.HI.f =0.005,
                     prev.W.f = 0.005,
                     prev.B.msf = 0.005,
                     prev.BI.msf =0.005,
                     prev.H.msf =0.005,
                     prev.HI.msf =0.005,
                     prev.W.msf = 0.5,
                     prev.B.msm = 0.1,
                     prev.BI.msm =0.1,
                     prev.H.msm =0.1,
                     prev.HI.msm =0.1,
                     prev.W.msm = 0.1,
                     prev.B.msmf = 0.05,
                     prev.BI.msmf =0.05,
                     prev.H.msmf =0.05,
                     prev.HI.msmf =0.05,
                     prev.W.msmf = 0.05,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  
  p$init.prev.age.slope.B.f <- prev.B.f / 38
  p$init.prev.age.slope.BI.f <- prev.BI.f / 38
  p$init.prev.age.slope.H.f <- prev.H.f / 38
  p$init.prev.age.slope.HI.f <- prev.HI.f / 38
  p$init.prev.age.slope.W.f <- prev.W.f / 38
  p$init.prev.age.slope.B.msf <- prev.B.msf / 38
  p$init.prev.age.slope.BI.msf <- prev.BI.msf / 38
  p$init.prev.age.slope.H.msf <- prev.H.msf / 38
  p$init.prev.age.slope.HI.msf <- prev.HI.msf / 38
  p$init.prev.age.slope.W.msf <- prev.W.msf / 38
  p$init.prev.age.slope.B.msm <- prev.B.msm / 38
  p$init.prev.age.slope.BI.msm <- prev.BI.msm / 38
  p$init.prev.age.slope.H.msm <- prev.H.msm / 38
  p$init.prev.age.slope.HI.msm <- prev.HI.msm / 38
  p$init.prev.age.slope.W.msm <- prev.W.msm / 38
  p$init.prev.age.slope.B.msmf <- prev.B.msmf / 38
  p$init.prev.age.slope.BI.msmf <- prev.BI.msmf / 38
  p$init.prev.age.slope.H.msmf <- prev.H.msmf / 38
  p$init.prev.age.slope.HI.msmf <- prev.HI.msmf / 38
  p$init.prev.age.slope.W.msmf <- prev.W.msmf / 38
 

  class(p) <- "init.net"
  return(p)
}


#' @title Epidemic Model Control Settings
#'
#' @description Sets the controls for stochastic network models simulated with
#'              \code{\link{netsim}}.
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes
#'        if used in conjunction with the \code{EpiModelHPC} package.
#' @param nsims Number of simulations.
#' @param ncores Number of cores per run, if parallelization is used within the
#'        \code{EpiModelHPC} package.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to 1 to run new
#'        simulation. This may also be set to 1 greater than the final time
#'        step of a previous simulation to resume the simulation with different
#'        parameters.
#' @param initialize.FUN Module function to use for initialization of the epidemic
#'        model.
#' @param aging.FUN Module function for aging.
#' @param deaths.FUN Module function for general and disease-realted deaths.
#' @param births.FUN Module function for births or entries into the population.
#' @param test.FUN Module function for diagnostic disease testing.
#' @param tx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param progress.FUN Module function for HIV disease progression.
#' @param vl.FUN Module function for HIV viral load evolution.
#' @param aiclass.FUN Module function for one-off AI risk class transitions.
#' @param roleclass.FUN Module function for transitions in sexual roles.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param disclose.FUN Module function for HIV status disclosure.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param riskhist.FUN Module function to calculate risk history for uninfected
#'        persons in the population.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param trans.FUN Module function to stochastically simulate disease transmission
#'        over acts given individual and dyadic attributes.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param verbose If \code{TRUE}, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param verbose.int Integer specifying the interval between time steps at which
#'        progress is printed.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{control_shamp}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
control_shamp <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        initialize.FUN = initialize_shamp,
                        aging.FUN = aging_msm,
                        deaths.FUN = deaths_shamp,
                        births.FUN = births_shamp,
                        test.FUN = test_shamp,
                        tx.FUN = tx_shamp,
                        prep.FUN = prep_shamp,
                        progress.FUN = progress_msm,
                        vl.FUN = vl_shamp,
                        #aiclass.FUN = NULL,
                        #roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_shamp,
                        disclose.FUN = disclose_shamp,
                        acts.FUN = acts_shamp,
                        condoms.FUN = condoms_shamp,
                        #riskhist.FUN = riskhist_shamp,
                        #position.FUN = position_shamp,
                        #trans.FUN = trans_msm,
                        prev.FUN = prevalence_shamp,
                        #verbose.FUN = verbose_msm,
                        save.nwstats = FALSE,
                        verbose = TRUE,
                        verbose.int = 1,
                        ...) {

  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)

  p$save.other = c("attr", "temp", "el", "p")

  p$save.network = FALSE

  class(p) <- "control.net"
  return(p)
}

