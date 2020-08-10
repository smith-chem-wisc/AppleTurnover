using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using System.Linq;
using System.Xaml;
using System.Diagnostics;
using FlashLFQ;
using MathNet.Numerics.Statistics;

namespace AppleTurnover
{
    public static class NonLinearRegression
    {
        //below values are arbitrary heuristics
        private const double ITERATIVE_SHIFT = 0.001; //original 0.00001
        private const double MAX_KST_VALUE = 2;
        private const double MAX_KBT_VALUE = 0.5;
        private const double MAX_KAO_VALUE = 5;

        private const double MIN_PARAMETER_VALUE = ITERATIVE_SHIFT * 10;
        //private const int NUM_TRAINING_GROUPS = 6;
        //private const int NUM_TRAINING_POINTS_PER_GROUP = 20;
        private const int NUM_TRAINING_POINTS = 100;

        public static PoolParameters RegressionAnalysis(List<PeptideTurnoverObject> peptides, string filePath, Settings settings)
        {
            string directory = Directory.GetParent(filePath).FullName;
            string filename = Path.GetFileNameWithoutExtension(filePath);

            //initial training set
            //List<PeptideTurnoverObject> peptidesToDetermineStartingParameters = new List<PeptideTurnoverObject>(); //peptides.GetRange(0, Math.Min(NUM_TRAINING_POINTS, peptides.Count));
            List<string> linesOfDifferentStarts = new List<string>();

            //starting values taken from figure 6b of Guan et al., 2011, Anal. Chem. "Compartment Modeling for Mammalian Protein Turnover Studies by Stable Isotope Metabolic Labeling"
            double ksto = 0.7;
            double kbto = 0.026;
            double kaoo = 2.0;

            //get approximate half-lives for each peptide
            UpdateKbi(ksto, kbto, kaoo, peptides);

            //split into 6 sections
            //kbi of >0.2, 0.2-0.1, 0.1-0.05, 0.05-0.025, 0.025-0.0125, <0.0125
            //half lives of <3.5, 3.5-6.9, 6.9-13.9, 13.9-27.7, 27.7-55.5, >55.5
            //for (int i = 0; i < NUM_TRAINING_GROUPS; i++)
            //{
            //    double currentMin = i == 0 ? 0 : 0.0125 * Math.Pow(2, i - 1);
            //    double currentMax = i == 5 ? double.PositiveInfinity : 0.0125 * Math.Pow(2, i);
            //    List<PeptideTurnoverObject> peptidesForThisSection = peptides.Where(x => x.Kbi <= currentMax && x.Kbi > currentMin).OrderByDescending(x => x.Timepoints.Length).ThenByDescending(x => x.TotalIntensity).ToList();
            //    int numPeptidesToAddFromThisSection = Math.Min(peptidesForThisSection.Count, NUM_TRAINING_POINTS_PER_GROUP);
            //    for (int j = 0; j < numPeptidesToAddFromThisSection; j++)
            //    {
            //        peptidesToDetermineStartingParameters.Add(peptidesForThisSection[j]);
            //    }
            //}

            //grab twice as many training points as desired
            List<PeptideTurnoverObject> innerQuartilePeptides = peptides.OrderBy(x => x.Kbi).ToList().GetRange(peptides.Count / 4, peptides.Count / 2).ToList();
            List<PeptideTurnoverObject> peptidesToDetermineStartingParameters = innerQuartilePeptides.OrderByDescending(x => x.Timepoints.Length).ThenByDescending(x => x.TotalIntensity).ToList();
            peptidesToDetermineStartingParameters = peptidesToDetermineStartingParameters.GetRange(0, Math.Min(NUM_TRAINING_POINTS, peptidesToDetermineStartingParameters.Count)); //grab a subset


            List<double> kstList = new List<double>();
            List<double> kbtList = new List<double>();
            List<double> kaoList = new List<double>();
            List<double> errors = new List<double>();

            //foreach training point
            foreach (PeptideTurnoverObject peptide in peptidesToDetermineStartingParameters)
            {
                PoolParameters parameters = new PoolParameters(ksto, kbto, kaoo);
                OptimizeFit(parameters, new List<PeptideTurnoverObject> { peptide });

                kstList.Add(parameters.Kst);
                kbtList.Add(parameters.Kbt);
                kaoList.Add(parameters.Kao);
                errors.Add(peptide.Error);
            }

            //Get median of the values
            PoolParameters bestVariables = new PoolParameters(kstList.Median(), kbtList.Median(), kaoList.Median());

            //train on test peptides
            OptimizeFit(bestVariables, peptidesToDetermineStartingParameters);

            //train on inner quartile peptides
            OptimizeFit(bestVariables, innerQuartilePeptides);

            double bestKst = bestVariables.Kst;
            double bestKbt = bestVariables.Kbt;
            double bestKao = bestVariables.Kao;

            //For each peptide, apply the kst, kbt, and kao, but optimize for the kbi
            UpdateKbi(bestKst, bestKbt, bestKao, peptides, ITERATIVE_SHIFT);
            //fine tune
            UpdateKbi(bestKst, bestKbt, bestKao, peptides, ITERATIVE_SHIFT / 10);
            UpdateKbi(bestKst, bestKbt, bestKao, peptides, ITERATIVE_SHIFT / 100);

            //remove messy peptides
            if (settings.RemoveMessyPeptides)
            {
                for (int i = peptides.Count - 1; i >= 0; i--)
                {
                    var peptide = peptides[i];
                    double[] predictions = PredictRelativeFractionUsingThreeCompartmentModel(bestKst, bestKbt, bestKao, peptide.Kbi, peptide.Timepoints);
                    double[] actualValues = peptide.RelativeFractions;
                    for (int j = 0; j < predictions.Length; j++)
                    {
                        if (Math.Abs(actualValues[j] - predictions[j]) > 0.1)
                        {
                            peptides.RemoveAt(i);
                            break;
                        }
                    }
                }
                //update inner quartile peptides
                innerQuartilePeptides = peptides.OrderBy(x => x.Kbi).ToList().GetRange(peptides.Count / 4, peptides.Count / 2).ToList();
            }

            //train on inner quartile peptides again
            OptimizeFit(bestVariables, innerQuartilePeptides);

            //grid analysis
            //create array of different ksts, kbts, and kaos to see if we're in a local minimum
            while (true)
            {
                double bestError = double.PositiveInfinity;
                double[] ratiosForIteration = new double[] { 0.1, 0.2, 0.33, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 2, 4 };
                for (int i = 0; i < ratiosForIteration.Length; i++)
                {
                    double kstCurrent = bestVariables.Kst * ratiosForIteration[i];
                    if (kstCurrent < MAX_KST_VALUE && kstCurrent > MIN_PARAMETER_VALUE)
                    {
                        for (int j = 0; j < ratiosForIteration.Length; j++)
                        {
                            double kbtCurrent = bestVariables.Kbt * ratiosForIteration[j];
                            if (kbtCurrent < MAX_KBT_VALUE && kbtCurrent > MIN_PARAMETER_VALUE)
                            {
                                for (int k = 0; k < ratiosForIteration.Length; k++)
                                {
                                    double kaoCurrent = bestVariables.Kao * ratiosForIteration[k];
                                    if (kaoCurrent < MAX_KAO_VALUE && kaoCurrent > MIN_PARAMETER_VALUE)
                                    {
                                        UpdateKbi(kstCurrent, kbtCurrent, kaoCurrent, innerQuartilePeptides, ITERATIVE_SHIFT);
                                        UpdateKbi(kstCurrent, kbtCurrent, kaoCurrent, innerQuartilePeptides, ITERATIVE_SHIFT / 10);
                                        UpdateKbi(kstCurrent, kbtCurrent, kaoCurrent, innerQuartilePeptides, ITERATIVE_SHIFT / 100);

                                        double currentError = innerQuartilePeptides.Sum(x => x.Error);
                                        if (currentError < bestError)
                                        {
                                            bestError = currentError;
                                            bestKst = kstCurrent;
                                            bestKbt = kbtCurrent;
                                            bestKao = kaoCurrent;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (bestVariables.Kst.Equals(bestKst) &&
                    bestVariables.Kbt.Equals(bestKbt) &&
                    bestVariables.Kao.Equals(bestKao))
                {
                    break;
                }
                else
                {
                    bestVariables.Kst = bestKst;
                    bestVariables.Kbt = bestKbt;
                    bestVariables.Kao = bestKao;

                    //train on all inner quartile peptides again
                    OptimizeFit(bestVariables, innerQuartilePeptides);
                }
            }

            bestKst = bestVariables.Kst;
            bestKbt = bestVariables.Kbt;
            bestKao = bestVariables.Kao;

            //For each peptide, apply the kst, kbt, and kao, and optimize for the kbi
            UpdateKbi(bestKst, bestKbt, bestKao, peptides, ITERATIVE_SHIFT);
            //fine tune
            UpdateKbi(bestKst, bestKbt, bestKao, peptides, ITERATIVE_SHIFT / 10);
            UpdateKbi(bestKst, bestKbt, bestKao, peptides, ITERATIVE_SHIFT / 100);

            //use the monte carlo method to estimate the 95% confidence interval
            Parallel.ForEach(Partitioner.Create(0, peptides.Count), new ParallelOptions { MaxDegreeOfParallelism = 8 },
    (partitionRange, loopState) =>
    {
        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
        {
            UpdateKbiConfidenceInterval(bestKst, bestKbt, bestKao, peptides[i], ITERATIVE_SHIFT);
        }
    });

            List<string> linesToWrite = new List<string>();
            linesToWrite.Add("Peptide\tProtein\tProteoform\tHalf-life\tLowerConfidenceInterval\tUpperConfidenceInterval\tError (MSE)\tSummed Intensity\tNumber of Ratios");
            string previousFullSeq = "";
            foreach (PeptideTurnoverObject peptide in peptides.OrderBy(x=>x.FullSequence))
            {
                if (!previousFullSeq.Equals(peptide.FullSequence))
                {
                    previousFullSeq = peptide.FullSequence;
                    linesToWrite.Add(peptide.FullSequence + '\t' + peptide.Protein + '\t'+peptide.Proteoform + '\t' + (Math.Log(2, Math.E) / peptide.Kbi).ToString() + '\t' + (Math.Log(2, Math.E) / peptide.HighKbi).ToString() + '\t' + (Math.Log(2, Math.E) / peptide.LowKbi).ToString() + '\t' + peptide.Error.ToString() + '\t' + peptide.TotalIntensity.ToString() + '\t' + peptide.Timepoints.Length.ToString());
                }
            }
            //output each peptide with its sequence, kbi, 95% confidence interval, and protein
            File.WriteAllLines(Path.Combine(directory, filename + "_PeptideTurnoverResults.tsv"), linesToWrite);


            return new PoolParameters(bestKst, bestKbt, bestKao);
        }

        public static void OptimizeFit(PoolParameters bestVariables, List<PeptideTurnoverObject> peptides)
        {
            bool optimizing = true;
            double updatedError = 0;
            bool increaseParameter = true;
            double kst = bestVariables.Kst;
            double kbt = bestVariables.Kbt;
            double kao = bestVariables.Kao;

            UpdateKbi(kst, kbt, kao, peptides, ITERATIVE_SHIFT);
            double previousError = peptides.Sum(x => x.Error);
            double bestError = previousError;

            while (optimizing)
            {
                optimizing = false;
                kst = bestVariables.Kst;
                kbt = bestVariables.Kbt;
                kao = bestVariables.Kao;
                PoolParameters variables = new PoolParameters(kst, kbt, kao);

                //KST
                if (kst + ITERATIVE_SHIFT < Math.Min(MAX_KST_VALUE, kao))
                {
                    UpdateKbi(kst + ITERATIVE_SHIFT, kbt, kao, peptides, ITERATIVE_SHIFT);
                    updatedError = peptides.Sum(x => x.Error);
                    increaseParameter = true;
                }
                else
                {
                    updatedError = double.PositiveInfinity;
                }
                if (!(updatedError < previousError) && kst - ITERATIVE_SHIFT > Math.Max(MIN_PARAMETER_VALUE, kbt))
                {
                    UpdateKbi(kst - ITERATIVE_SHIFT, kbt, kao, peptides, ITERATIVE_SHIFT);
                    updatedError = peptides.Sum(x => x.Error);
                    increaseParameter = false;
                }

                if (previousError > updatedError)
                {
                    optimizing = true;
                    double diff = (previousError - updatedError) / ITERATIVE_SHIFT;
                    if (diff > ITERATIVE_SHIFT * 3)
                    {
                        diff = Math.Round(diff / ITERATIVE_SHIFT) * ITERATIVE_SHIFT;
                        double tempError = double.PositiveInfinity;
                        if (increaseParameter)
                        {
                            if (kst + diff < Math.Min(MAX_KST_VALUE, kao)) //kst should be less than kao
                            {
                                UpdateKbi(kst + diff, kbt, kao, peptides, ITERATIVE_SHIFT);
                                tempError = peptides.Sum(x => x.Error);
                            }
                        }
                        else
                        {
                            if (kst - diff > Math.Max(MIN_PARAMETER_VALUE, kbt))
                            {
                                UpdateKbi(kst - diff, kbt, kao, peptides, ITERATIVE_SHIFT);
                                tempError = peptides.Sum(x => x.Error);
                            }
                        }

                        if (tempError > updatedError)
                        {
                            diff = ITERATIVE_SHIFT;
                        }
                        else
                        {
                            updatedError = tempError;
                        }
                    }
                    else
                    {
                        diff = ITERATIVE_SHIFT;
                    }

                    if (increaseParameter)
                    {
                        variables.Kst += diff;
                    }
                    else
                    {
                        variables.Kst -= diff;
                    }

                    if (updatedError < bestError)
                    {
                        bestError = updatedError;
                        bestVariables.Kst = variables.Kst;
                        bestVariables.Kbt = kbt;
                        bestVariables.Kao = kao;
                    }
                }

                //KBT
                if (kbt + ITERATIVE_SHIFT < Math.Min(MAX_KBT_VALUE, kst))
                {
                    UpdateKbi(kst, kbt + ITERATIVE_SHIFT, kao, peptides, ITERATIVE_SHIFT);
                    updatedError = peptides.Sum(x => x.Error);
                    increaseParameter = true;
                }
                else
                {
                    updatedError = double.PositiveInfinity;
                }
                if (!(updatedError < previousError) && kbt - ITERATIVE_SHIFT > MIN_PARAMETER_VALUE)
                {
                    UpdateKbi(kst, kbt - ITERATIVE_SHIFT, kao, peptides, ITERATIVE_SHIFT);
                    updatedError = peptides.Sum(x => x.Error);
                    increaseParameter = false;
                }

                if (previousError > updatedError)
                {
                    optimizing = true;
                    double diff = (previousError - updatedError) / ITERATIVE_SHIFT;
                    if (diff > ITERATIVE_SHIFT * 3)
                    {
                        diff = Math.Round(diff / ITERATIVE_SHIFT) * ITERATIVE_SHIFT;
                        double tempError = double.PositiveInfinity;
                        if (increaseParameter)
                        {
                            if (kbt + diff < Math.Min(MAX_KBT_VALUE, kst)) //kbt should be less than kst
                            {
                                UpdateKbi(kst, kbt + diff, kao, peptides, ITERATIVE_SHIFT);
                                tempError = peptides.Sum(x => x.Error);
                            }
                        }
                        else
                        {
                            if (kbt - diff > MIN_PARAMETER_VALUE)
                            {
                                UpdateKbi(kst, kbt - diff, kao, peptides, ITERATIVE_SHIFT);
                                tempError = peptides.Sum(x => x.Error);
                            }
                        }

                        if (tempError > updatedError)
                        {
                            diff = ITERATIVE_SHIFT;
                        }
                        else
                        {
                            updatedError = tempError;
                        }
                    }
                    else
                    {
                        diff = ITERATIVE_SHIFT;
                    }

                    if (increaseParameter)
                    {
                        variables.Kbt += diff;
                    }
                    else
                    {
                        variables.Kbt -= diff;
                    }

                    if (updatedError < bestError)
                    {
                        bestError = updatedError;
                        bestVariables.Kst = kst;
                        bestVariables.Kbt = variables.Kbt;
                        bestVariables.Kao = kao;
                    }
                }

                //KAO
                if (kao + ITERATIVE_SHIFT < MAX_KAO_VALUE)
                {
                    UpdateKbi(kst, kbt, kao + ITERATIVE_SHIFT, peptides, ITERATIVE_SHIFT);
                    updatedError = peptides.Sum(x => x.Error);
                    increaseParameter = true;
                }
                else
                {
                    updatedError = double.PositiveInfinity;
                }
                if (!(updatedError < previousError) && kao - ITERATIVE_SHIFT > Math.Max(MIN_PARAMETER_VALUE, kst)) //kao should be greater than kst
                {
                    UpdateKbi(kst, kbt, kao - ITERATIVE_SHIFT, peptides, ITERATIVE_SHIFT);
                    updatedError = peptides.Sum(x => x.Error);
                    increaseParameter = false;
                }

                if (previousError > updatedError)
                {
                    optimizing = true;
                    double diff = (previousError - updatedError) / ITERATIVE_SHIFT;
                    if (diff > ITERATIVE_SHIFT * 3)
                    {
                        diff = Math.Round(diff / ITERATIVE_SHIFT) * ITERATIVE_SHIFT;
                        double tempError = double.PositiveInfinity;
                        if (increaseParameter)
                        {
                            if (kao + diff < MAX_KAO_VALUE)
                            {
                                UpdateKbi(kst, kbt, kao + diff, peptides, ITERATIVE_SHIFT);
                                tempError = peptides.Sum(x => x.Error);
                            }
                        }
                        else
                        {
                            if (kao - diff > Math.Max(MIN_PARAMETER_VALUE, kst)) //kao should be greater than kst
                            {
                                UpdateKbi(kst, kbt, kao - diff, peptides, ITERATIVE_SHIFT);
                                tempError = peptides.Sum(x => x.Error);
                            }
                        }

                        if (tempError > updatedError)
                        {
                            diff = ITERATIVE_SHIFT;
                        }
                        else
                        {
                            updatedError = tempError;
                        }
                    }
                    else
                    {
                        diff = ITERATIVE_SHIFT;
                    }

                    if (increaseParameter)
                    {
                        variables.Kao += diff;
                    }
                    else
                    {
                        variables.Kao -= diff;
                    }

                    if (updatedError < bestError)
                    {
                        bestError = updatedError;
                        bestVariables.Kst = kst;
                        bestVariables.Kbt = kbt;
                        bestVariables.Kao = variables.Kao;
                    }
                }

                //test the new gradient to ensure we don't get trapped
                UpdateKbi(variables.Kst, variables.Kbt, variables.Kao, peptides, ITERATIVE_SHIFT);
                updatedError = peptides.Sum(x => x.Error);
                if (updatedError < bestError)
                {
                    bestError = updatedError;
                    bestVariables.Kst = variables.Kst;
                    bestVariables.Kbt = variables.Kbt;
                    bestVariables.Kao = variables.Kao;
                }
                else //reset to the old
                {
                    UpdateKbi(bestVariables.Kst, bestVariables.Kbt, bestVariables.Kao, peptides, ITERATIVE_SHIFT);
                }

                previousError = bestError;
            }
        }

        public static void UpdateKbi(double kst, double kbt, double kao, List<PeptideTurnoverObject> peptides, double ITERATIVE_SHIFT = ITERATIVE_SHIFT)
        {
            if (peptides.Count < 100) //if not worth parallelizing
            {
                for (int i = 0; i < peptides.Count; i++)
                {
                    UpdateKbi(kst, kbt, kao, peptides[i], ITERATIVE_SHIFT);
                }
            }
            else
            {
                Parallel.ForEach(Partitioner.Create(0, peptides.Count), new ParallelOptions { MaxDegreeOfParallelism = 8 },
                    (partitionRange, loopState) =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            UpdateKbi(kst, kbt, kao, peptides[i], ITERATIVE_SHIFT);
                        }
                    });
            }
        }

        public static void UpdateKbi(double kst, double kbt, double kao, PeptideTurnoverObject peptide, double ITERATIVE_SHIFT = ITERATIVE_SHIFT)
        {
            //List<PeptideTurnoverValues> peptidesToRemove = new List<PeptideTurnoverValues>();
            const double MAX_KBI = 2.0;
            double MIN_KBI = ITERATIVE_SHIFT * 10;
            double ARBITRARY_GRADIENT_FACTOR = 50 * ITERATIVE_SHIFT;

            double[] timepoints = peptide.Timepoints;
            double[] relativeFractions = peptide.RelativeFractions;
            double kbi = peptide.Kbi;
            // if (peptide.Error==double.PositiveInfinity) //first time
            {
                peptide.Error = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi, timepoints, relativeFractions);
                peptide.TemporaryError = peptide.Error;
            }
            double originalError = peptide.Error;

            double updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi + ITERATIVE_SHIFT, timepoints, relativeFractions);
            bool increaseKbi = true;
            if (!(updatedError < originalError))
            {
                updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi - ITERATIVE_SHIFT, timepoints, relativeFractions);
                increaseKbi = false;
            }

            if (originalError > updatedError)
            {
                double diff = (originalError - updatedError) / ARBITRARY_GRADIENT_FACTOR;
                if (diff > ITERATIVE_SHIFT * 3)
                {
                    diff = Math.Round(diff / ITERATIVE_SHIFT) * ITERATIVE_SHIFT;
                    double tempError;
                    if (increaseKbi)
                    {
                        if (kbi + diff > MAX_KBI)
                        {
                            diff = MAX_KBI - diff - ITERATIVE_SHIFT;
                        }
                        tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi + diff, timepoints, relativeFractions);
                    }
                    else
                    {
                        if (kbi - diff < MIN_KBI)
                        {
                            diff = kbi - MIN_KBI - ITERATIVE_SHIFT;
                        }
                        tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi - diff, timepoints, relativeFractions);
                    }

                    if (!(tempError < updatedError))
                    {
                        diff = ITERATIVE_SHIFT;
                    }
                    else
                    {
                        updatedError = tempError;
                    }
                }
                else
                {
                    diff = ITERATIVE_SHIFT;
                }

                if (increaseKbi)
                {
                    kbi += diff;
                }
                else
                {
                    kbi -= diff;
                }

                while (originalError > updatedError)
                {
                    if (kbi > MAX_KBI || kbi < MIN_KBI) //max and min allowed
                    {
                        //peptidesToRemove.Add(peptide);
                        break;
                    }
                    originalError = updatedError;

                    updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi + ITERATIVE_SHIFT, timepoints, relativeFractions);
                    increaseKbi = true;
                    if (!(updatedError < originalError))
                    {
                        updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi - ITERATIVE_SHIFT, timepoints, relativeFractions);
                        increaseKbi = false;
                    }

                    diff = (originalError - updatedError) / ARBITRARY_GRADIENT_FACTOR;
                    if (diff > ITERATIVE_SHIFT * 3)
                    {
                        diff = Math.Round(diff / ITERATIVE_SHIFT) * ITERATIVE_SHIFT;
                        double tempError;
                        if (increaseKbi)
                        {
                            if (kbi + diff > MAX_KBI)
                            {
                                diff = MAX_KBI - diff - ITERATIVE_SHIFT;
                            }
                            tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi + diff, timepoints, relativeFractions);
                        }
                        else
                        {
                            if (kbi - diff < MIN_KBI)
                            {
                                diff = kbi - MIN_KBI - ITERATIVE_SHIFT;
                            }
                            tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, kao, kbi - diff, timepoints, relativeFractions);
                        }

                        if (!(tempError < updatedError))
                        {
                            diff = ITERATIVE_SHIFT;
                        }
                        else
                        {
                            updatedError = tempError;
                        }
                    }
                    else
                    {
                        diff = ITERATIVE_SHIFT;
                    }

                    if (increaseKbi)
                    {
                        kbi += diff;
                    }
                    else
                    {
                        kbi -= diff;
                    }
                }
                peptide.TemporaryError = originalError;

                if (increaseKbi)
                {
                    kbi -= diff;
                }
                else
                {
                    kbi += diff;
                }
                peptide.Kbi = kbi;
            }

            peptide.UpdateError();

        }

        /// <param name="kst">Constant for amino acid compartment to protein compartment</param>
        /// <param name="kbt">Constant for protein compartment to amino acid compartment</param>
        /// <param name="kao">Constant for amino acid compartment to waste</param>
        /// <param name="kbi">Constant for the protein of interest's degredation</param>
        private static double CalculateErrorForThreeCompartmentModelFit(double kst, double kbt, double kao, double kbi, double[] timepoints, double[] relativeFraction)
        {
            double[] residuals = new double[timepoints.Length];
            double[] predictedValues = PredictRelativeFractionUsingThreeCompartmentModel(kst, kbt, kao, kbi, timepoints);

            for (int i = 0; i < timepoints.Length; i++)
            {
                residuals[i] = relativeFraction[i] - predictedValues[i];
            }

            //SSE is the sum, this is the average
            return residuals.Average(x => Math.Pow(x, 2));
        }

        //assumes timepoints are sorted
        public static double[] PredictRelativeFractionUsingThreeCompartmentModel(double kst, double kbt, double kao, double kbi, double[] timepoints)
        {
            double sumOfStOaBt = kst + kbt + kao;
            double sqrtForUV = Math.Sqrt(Math.Pow(sumOfStOaBt, 2) - 4 * kao * kbt);
            double u = (sumOfStOaBt - sqrtForUV) / 2;
            double v = (sumOfStOaBt + sqrtForUV) / 2;
            double yu = kao * kbi * (u - kbt) / ((u - v) * (u - kbi) * u);
            double yv = kao * kbi * (v - kbt) / ((v - u) * (v - kbi) * v);
            double ykbi = kao * (kbi - kbt) / ((u - kbi) * (v - kbi));

            double[] predictedValues = new double[timepoints.Length];
            double previousTimepoint = double.NegativeInfinity;
            for (int i = 0; i < timepoints.Length; i++)
            {
                double t = Math.Max(timepoints[i], 0);
                if (!t.Equals(previousTimepoint))
                {
                    previousTimepoint = t;
                    predictedValues[i] = 1 + yu * Math.Exp(-u * t) + yv * Math.Exp(-v * t) + ykbi * Math.Exp(-kbi * t);
                }
                else
                {
                    predictedValues[i] = predictedValues[i - 1];
                }
            }
            return predictedValues;
        }

        public static void UpdateKbiConfidenceInterval(double kst, double kbt, double kao, PeptideTurnoverObject peptide, double ITERATIVE_SHIFT)
        {
            double[] timepoints = peptide.Timepoints;
            //int numMeasurements = timepoints.Length;
            double[] relativeFraction = peptide.RelativeFractions;
            //double[] residuals = new double[numMeasurements];
            double originalKbi = peptide.Kbi; //save it so it doesn't get overwritten
            //double[] distinctPredictedValues = PredictRelativeFractionUsingThreeCompartmentModel(kst, kbt, kao, peptide.Kbi, timepoints).Distinct().ToArray();
            //double[] uniqueTimePoints = timepoints.Distinct().ToArray();
            List<List<double>> relativeFractionsForEachTimePoint = new List<List<double>>();
            double previousTimepoint = timepoints[0];

            int previousI = 0;
            //get predicted values
            double[] predicted = PredictRelativeFractionUsingThreeCompartmentModel(kst, kbt, kao, peptide.Kbi, timepoints);
            //find sigma for each timepoint

            List<double> sigmasForEachTimePoint = new List<double>();
            List<double> sigmasForThisTimepoint = new List<double> { Math.Pow((relativeFraction[0] - predicted[0]), 2) };
            List<double> relativeFractionsForThisTimepoint = new List<double> { relativeFraction[0] };
            for (int i = 1; i < timepoints.Length; i++)
            {
                if (timepoints[i].Equals(previousTimepoint))
                {
                    sigmasForThisTimepoint.Add(Math.Pow((relativeFraction[i] - predicted[i]), 2));
                    relativeFractionsForThisTimepoint.Add(relativeFraction[i]);
                }
                else
                {
                    for (; previousI < i; previousI++)
                    {
                        sigmasForEachTimePoint.Add(Math.Sqrt(sigmasForThisTimepoint.Sum() / (Math.Max(sigmasForThisTimepoint.Count - 1, 1)))); //prevent infinity if only one point
                        relativeFractionsForEachTimePoint.Add(relativeFractionsForThisTimepoint);
                    }
                    sigmasForThisTimepoint = new List<double> { Math.Pow((relativeFraction[i] - predicted[i]), 2) };
                    relativeFractionsForThisTimepoint = new List<double> { relativeFraction[i] };
                    previousTimepoint = timepoints[i];
                }
            }
            for (; previousI < timepoints.Length; previousI++)
            {
                sigmasForEachTimePoint.Add(Math.Sqrt(sigmasForThisTimepoint.Sum() / (Math.Max(sigmasForThisTimepoint.Count - 1, 1)))); //prevent infinity if only one point
                relativeFractionsForEachTimePoint.Add(relativeFractionsForThisTimepoint); //add the last one
            }
            ////BELOW CODE ONLY CREATES ONE LIST FOR EACH UNIQUE TIMEPOINT
            //List<double> relativeFractionsForThisTimepoint = new List<double> { relativeFraction[0] };
            //for(int i=1; i<timepoints.Length;i++)
            //{
            //    if(timepoints[i].Equals(previousTimepoint))
            //    {
            //        relativeFractionsForThisTimepoint.Add(relativeFraction[i]);
            //    }
            //    else
            //    {
            //        relativeFractionsForEachTimePoint.Add(relativeFractionsForThisTimepoint);
            //        relativeFractionsForThisTimepoint = new List<double> { relativeFraction[i] };
            //        previousTimepoint = timepoints[i];
            //    }
            //}
            //relativeFractionsForEachTimePoint.Add(relativeFractionsForThisTimepoint); //add the last one

            //use as an index
            double originalError = peptide.Error;
            //for (int i = 0; i < numMeasurements; i++)
            //{
            //    residuals[i] = relativeFraction[i] - predictedValues[i];
            //}

            int NUM_SIMULATIONS = 200;
            int lowPercentile = (int)Math.Round((NUM_SIMULATIONS-1) * 2.5 / 100); //2.5 percentile
            int highPercentile = (int)Math.Round((NUM_SIMULATIONS-1) * 97.5 / 100); //97.5 percentile

            //simulate data through monte carlo (or bootstrapping, commented out)
            Random rng = new Random(1);
            //double sigma = Math.Sqrt(peptide.Error * timepoints.Length / (timepoints.Length - 4)); //The error is actually the MSE, we want SSE, so multiply, by timepoints, divide by degrees of freedom with 4 free variables
            double[] bootstrapKbis = new double[NUM_SIMULATIONS];
            for (int s = 0; s < NUM_SIMULATIONS; s++)
            {
                double[] simulatedData = new double[timepoints.Length];
                for (int i = 0; i < timepoints.Length; i++)
                {
                    //bootstrap
                    //simulatedData[i] = predictedValues[i] + residuals[rng.Next(numMeasurements)];
                    //monte carlo
                    List<double> actualValuesForThisTimepoint = relativeFractionsForEachTimePoint[i];
                    //take an actual value (sampling with replacement) and then add variance based on the residual
                    double sampleValue = actualValuesForThisTimepoint[rng.Next(actualValuesForThisTimepoint.Count)];
                    double varianceToAdd = NormInv(rng.NextDouble(), 0, sigmasForEachTimePoint[i]);
                    double simulatedValue = sampleValue + varianceToAdd;
                    //don't allow the simulated value to be an unreal value
                    if (simulatedValue > 1)
                    {
                        simulatedValue = 1;
                    }
                    if (simulatedValue < 0)
                    {
                        simulatedValue = 0;
                    }
                    simulatedData[i] = simulatedValue;
                    //simulatedData[i] = distinctPredictedValues[i] + NormInv(rng.NextDouble(), 0, sigma);
                }
                peptide.RelativeFractions = simulatedData;

                //optimize on the simulated data
                UpdateKbi(kst, kbt, kao, peptide, ITERATIVE_SHIFT);
                UpdateKbi(kst, kbt, kao, peptide, ITERATIVE_SHIFT / 10);
                UpdateKbi(kst, kbt, kao, peptide, ITERATIVE_SHIFT / 100);
                bootstrapKbis[s] = peptide.Kbi;
            }
            ////THE BELOW COMMENTED CODE ALLOWED FOR ONLY ONE SIMULATED POINT PER TIMEPOINT
            //double sigma = Math.Sqrt(peptide.Error * timepoints.Length / (timepoints.Length - 4)); //The error is actually the MSE, we want SSE, so multiply, by timepoints, divide by degrees of freedom with 4 free variables
            //double[] bootstrapKbis = new double[NUM_SIMULATIONS];
            //peptide.Timepoints = uniqueTimePoints;
            //for (int s = 0; s < NUM_SIMULATIONS; s++)
            //{
            //    double[] simulatedData = new double[uniqueTimePoints.Length];
            //    for (int i = 0; i < uniqueTimePoints.Length; i++)
            //    {
            //        //bootstrap
            //        //simulatedData[i] = predictedValues[i] + residuals[rng.Next(numMeasurements)];
            //        //monte carlo
            //        List<double> actualValuesForThisTimepoint = relativeFractionsForEachTimePoint[i];
            //        simulatedData[i] = actualValuesForThisTimepoint[rng.Next(actualValuesForThisTimepoint.Count)] + NormInv(rng.NextDouble(), 0, sigma);
            //        //simulatedData[i] = distinctPredictedValues[i] + NormInv(rng.NextDouble(), 0, sigma);
            //    }
            //    peptide.RelativeFractions = simulatedData;

            //    //optimize on the simulated data
            //    UpdateKbi(kst, kbt, kao, peptide, ITERATIVE_SHIFT);
            //    bootstrapKbis[s] = peptide.Kbi;
            //}
            peptide.Kbi = originalKbi;
            peptide.RelativeFractions = relativeFraction;
            peptide.Timepoints = timepoints;
            peptide.Error = originalError;
            Array.Sort(bootstrapKbis);
            peptide.LowKbi = bootstrapKbis[lowPercentile]; //2.5 percentile
            peptide.HighKbi = bootstrapKbis[highPercentile]; //97.5 percentile
            peptide.MonteCarloKbis = bootstrapKbis;
        }

        public static List<PeptideTurnoverObject> GetProteinInfo(List<PeptideTurnoverObject> peptides, string filePath, List<IGrouping<string, PeptideTurnoverObject>> grouping, string analysisType)
        {
            //group peptides by protein
            List<PeptideTurnoverObject> proteinsToReturn = new List<PeptideTurnoverObject>();
            foreach (var group in grouping)
            {
                string proteinName = group.Key;
                var peptidesForThisProtein = group.OrderBy(x => x.BaseSequence).ThenBy(x => x.FullSequence).ToList();
                List<double> kbis = new List<double>();
                List<double> timepoints = new List<double>();
                List<double> rfs = new List<double>();
                List<string> filenames = new List<string>();
                List<double> intensities = new List<double>();
                double intensity = 0;
                string allPeptideSequences = "";
                foreach (var peptide in peptidesForThisProtein)
                {
                    kbis.AddRange(peptide.MonteCarloKbis);
                    timepoints.AddRange(peptide.Timepoints);
                    rfs.AddRange(peptide.RelativeFractions);
                    filenames.AddRange(peptide.Filenames);
                    intensities.AddRange(peptide.Intensities);
                    intensity += peptide.TotalIntensity;
                    allPeptideSequences += peptide.FullSequence + ";";
                }
                //keep excel compatible
                if (allPeptideSequences.Length > 10000)
                {
                    allPeptideSequences = "Too Many Sequences";
                }
                kbis.Sort();
                //get important measurements
                PeptideTurnoverObject protein = new PeptideTurnoverObject(allPeptideSequences, timepoints.ToArray(), rfs.ToArray(),
                    filenames.ToArray(), intensities.ToArray(), intensity, peptides.First().FileName, proteinName);

                //we do 200*n iterations, so we want the average half-life (not kbi) of the 99*nth and 100*nth index to get the median half-life. Then transform back to kbi
                double medianHalfLifeLow = Math.Log(2, Math.E) / kbis[kbis.Count / 2 - 1];
                double medianHalfLifeHigh = Math.Log(2, Math.E) / kbis[kbis.Count / 2];
                double medianHalfLife = (medianHalfLifeHigh + medianHalfLifeLow) / 2;
                protein.Kbi = Math.Log(2, Math.E) / medianHalfLife; 
                protein.LowKbi = kbis[(int)Math.Round((kbis.Count-1) * 2.5 / 100)]; //2.5 percentile
                protein.HighKbi = kbis[(int)Math.Round((kbis.Count - 1) * 97.5 / 100)]; //97.5 percentile
                protein.MonteCarloKbis = kbis.ToArray();
                protein.Error = peptidesForThisProtein.Count; //not the actual error
                proteinsToReturn.Add(protein);
            }

            string directory = Directory.GetParent(filePath).FullName;
            string filename = Path.GetFileNameWithoutExtension(filePath);

            List<string> linesToWrite = new List<string>();
            linesToWrite.Add("Protein\tHalf Life\tLowerConfidenceInterval\tUpperConfidenceInterval\tSummed Intensity\tNumber of Ratios\tNumber of Peptides\tPeptideSequences");
            foreach (PeptideTurnoverObject protein in proteinsToReturn)
            {
                linesToWrite.Add(protein.Protein + '\t' + (Math.Log(2, Math.E) / protein.Kbi).ToString() + '\t' +
                    (Math.Log(2, Math.E) / protein.HighKbi).ToString() + '\t' +
                    (Math.Log(2, Math.E) / protein.LowKbi).ToString() + '\t' + protein.TotalIntensity.ToString() + '\t' +
                    protein.Timepoints.Length.ToString() + '\t' + protein.Error.ToString() + '\t' + protein.FullSequence);
            }

            //output each peptide with its sequence, kbi, 95% confidence interval, and protein
            File.WriteAllLines(Path.Combine(directory, filename + "_" + analysisType + "TurnoverResults.tsv"), linesToWrite);

            return proteinsToReturn;
        }

        /// <summary>
        /// Given a probability, a mean, and a standard deviation, an x value can be calculated.
        /// https://stackoverflow.com/questions/2901750/is-there-a-c-sharp-library-that-will-perform-the-excel-norminv-function
        /// </summary>
        /// <returns></returns>
        public static double NormInv(double probability)
        {
            const double a1 = -39.6968302866538;
            const double a2 = 220.946098424521;
            const double a3 = -275.928510446969;
            const double a4 = 138.357751867269;
            const double a5 = -30.6647980661472;
            const double a6 = 2.50662827745924;

            const double b1 = -54.4760987982241;
            const double b2 = 161.585836858041;
            const double b3 = -155.698979859887;
            const double b4 = 66.8013118877197;
            const double b5 = -13.2806815528857;

            const double c1 = -7.78489400243029E-03;
            const double c2 = -0.322396458041136;
            const double c3 = -2.40075827716184;
            const double c4 = -2.54973253934373;
            const double c5 = 4.37466414146497;
            const double c6 = 2.93816398269878;

            const double d1 = 7.78469570904146E-03;
            const double d2 = 0.32246712907004;
            const double d3 = 2.445134137143;
            const double d4 = 3.75440866190742;

            //Define break-points
            const double pLow = 0.02425;

            const double pHigh = 1 - pLow;

            //Define work variables
            double q;
            double result = 0;

            // if argument out of bounds.
            // set it to a value within desired precision.
            if (probability <= 0)
                probability = pLow;

            if (probability >= 1)
                probability = pHigh;

            if (probability < pLow)
            {
                //Rational approximation for lower region
                q = Math.Sqrt(-2 * Math.Log(probability));
                result = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
            }
            else if (probability <= pHigh)
            {
                //Rational approximation for lower region
                q = probability - 0.5;
                double r = q * q;
                result = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
                         (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
            }
            else if (probability < 1)
            {
                //Rational approximation for upper region
                q = Math.Sqrt(-2 * Math.Log(1 - probability));
                result = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
            }

            return result;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="probability"></param>
        /// <param name="mean"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double NormInv(double probability, double mean, double sigma)
        {
            double x = NormInv(probability);
            return sigma * x + mean;
        }
    }
}