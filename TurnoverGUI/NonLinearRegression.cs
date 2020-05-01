using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using System.Linq;

namespace TurnoverGUI
{
    public class NonLinearRegression
    {
        private static readonly string ParameterSaveFile = "_FittedParameters.txt";


        //below values are arbitrary heuristics
        private const double ITERATIVE_SHIFT = 0.001; //original 0.00001
        private const double MAX_KST_VALUE = 2;
        private const double MAX_KBT_VALUE = 0.5;
        private const double MAX_KOA_VALUE = 5;
        private const double MIN_PARAMETER_VALUE = ITERATIVE_SHIFT * 10;
        private const int NUM_TRAINING_POINTS = 100;

        public static PoolParameters RegressionAnalysis(List<PeptideTurnoverObject> peptides, string filePath)
        {
            //check if this has already been run
            string directory = Directory.GetParent(filePath).FullName;
            string filename = Path.GetFileNameWithoutExtension(filePath);
            string pathForFittedParameters = Path.Combine(directory, filename + ParameterSaveFile);

            double bestKst = 0;
            double bestKbt = 0;
            double bestKoa = 0;
            if (File.Exists(pathForFittedParameters))
            {
                string[] paramLine = File.ReadAllLines(pathForFittedParameters)[1].Split('\t').ToArray();
                bestKst = Convert.ToDouble(paramLine[0]);
                bestKbt = Convert.ToDouble(paramLine[1]);
                bestKoa = Convert.ToDouble(paramLine[2]);
            }
            else
            {
                //initial training set
                List<PeptideTurnoverObject> trainingPeptides = peptides.GetRange(0, Math.Min(NUM_TRAINING_POINTS, peptides.Count));

                List<string> linesOfDifferentStarts = new List<string>();
                //starting values taken from figure 6b of Guan et al., 2011, Anal. Chem. "Compartment Modeling for Mammalian Protein Turnover Studies by Stable Isotope Metabolic Labeling"
                double ksto = 0.7;
                double kbto = 0.026;
                double koao = 2.0;
                List<double> kstList = new List<double>();
                List<double> kbtList = new List<double>();
                List<double> koaList = new List<double>();
                List<double> errors = new List<double>();

                //foreach training point
                for (int index = 0; index < trainingPeptides.Count; index++)
                {
                    PeptideTurnoverObject peptide = trainingPeptides[index];
                    double[] bestVariables = new double[] { ksto, kbto, koao };
                    OptimizeFit(bestVariables, ksto, kbto, koao, peptide);

                    kstList.Add(bestVariables[0]);
                    kbtList.Add(bestVariables[1]);
                    koaList.Add(bestVariables[2]);
                    errors.Add(peptide.Error);
                }

                //Remove half of the training set that had the lowest errors.
                //Add additional peptides (required to have error under the median of the initial set) until we have 100 training points
                errors.Sort();
                double medianError = errors[errors.Count / 2];
                for (int i = NUM_TRAINING_POINTS - 1; i >= 0; i--)
                {
                    if (trainingPeptides[i].Error > medianError)
                    {
                        trainingPeptides.RemoveAt(i);
                        kstList.RemoveAt(i);
                        kbtList.RemoveAt(i);
                        koaList.RemoveAt(i);
                    }
                }

                for (int index = NUM_TRAINING_POINTS; index < peptides.Count; index++)
                {
                    if (trainingPeptides.Count == NUM_TRAINING_POINTS)
                    {
                        break;
                    }
                    PeptideTurnoverObject peptide = peptides[index];
                    double[] bestVariables = new double[] { ksto, kbto, koao };
                    OptimizeFit(bestVariables, ksto, kbto, koao, peptide);

                    //save the fits
                    if (peptide.Error < medianError)
                    {
                        trainingPeptides.Add(peptide);
                        kstList.Add(bestVariables[0]);
                        kbtList.Add(bestVariables[1]);
                        koaList.Add(bestVariables[2]);
                    }
                }


                //iterate through each of the models with all of the training peptides to determine the best error
                int bestIndex = 0;
                double bestModelError = double.PositiveInfinity;
                for (int modelIndex = 0; modelIndex < trainingPeptides.Count; modelIndex++)
                {
                    double kst = kstList[modelIndex];
                    double kbt = kbtList[modelIndex];
                    double koa = koaList[modelIndex];
                    if (kst < MAX_KST_VALUE - ITERATIVE_SHIFT && kbt < MAX_KBT_VALUE - ITERATIVE_SHIFT && koa < MAX_KOA_VALUE - ITERATIVE_SHIFT)
                    {
                        UpdateKbi(kst, kbt, koa, trainingPeptides, ITERATIVE_SHIFT);
                        double error = trainingPeptides.Sum(x => x.Error);
                        if (error < bestModelError)
                        {
                            bestModelError = error;
                            bestIndex = modelIndex;
                        }
                    }
                }

                //train all peptides simultaneously using the current best fit as the starting fit
                List<double> shiftsToUse = new List<double>();
                for (int i = 0; i < 1; i++)
                {
                    shiftsToUse.Add(ITERATIVE_SHIFT / (Math.Pow(10, i)));
                }
                foreach (double shift in shiftsToUse)
                {
                    double kst = kstList[bestIndex];
                    double kbt = kbtList[bestIndex];
                    double koa = koaList[bestIndex];
                    double[] variables = new double[] { kst, kbt, koa };

                    bool optimizing = true;
                    List<string> gradientMove = new List<string>();
                    double previousError = 0;
                    double updatedError = 0;
                    bool increaseParameter = true;
                    UpdateKbi(kst, kbt, koa, trainingPeptides, shift);
                    previousError = trainingPeptides.Sum(x => x.Error);
                    double bestError = previousError;
                    double[] bestVariables = new double[] { kst, kbt, koa };

                    //Optimizing fits across multiple peptides, not just a single one
                    while (optimizing)
                    {
                        optimizing = false;
                        kst = bestVariables[0];
                        kbt = bestVariables[1];
                        koa = bestVariables[2];
                        variables = new double[] { kst, kbt, koa };

                        //gradientMove.Add(kst.ToString() + '\t' + kbt.ToString() + '\t' + koa.ToString() + '\t' + previousError.ToString());
                        //if(gradientMove.Count>10000)
                        //{ }
                        //KST
                        if (kst + shift < Math.Min(MAX_KST_VALUE, koa))
                        {
                            UpdateKbi(kst + shift, kbt, koa, trainingPeptides, shift);
                            updatedError = trainingPeptides.Sum(x => x.Error);
                            increaseParameter = true;
                        }
                        else
                        {
                            updatedError = double.PositiveInfinity;
                        }
                        if (!(updatedError < previousError) && kst - shift > Math.Max(MIN_PARAMETER_VALUE, kbt))
                        {
                            UpdateKbi(kst - shift, kbt, koa, trainingPeptides, shift);
                            updatedError = trainingPeptides.Sum(x => x.Error);
                            increaseParameter = false;
                        }

                        if (previousError > updatedError)
                        {
                            optimizing = true;
                            double diff = (previousError - updatedError) / shift;
                            if (diff > shift * 3)
                            {
                                diff = Math.Round(diff / shift) * shift;
                                double tempError = double.PositiveInfinity;
                                if (increaseParameter)
                                {
                                    if (kst + diff < Math.Min(MAX_KST_VALUE, koa))
                                    {
                                        UpdateKbi(kst + diff, kbt, koa, trainingPeptides, shift);
                                        tempError = trainingPeptides.Sum(x => x.Error);
                                    }
                                }
                                else
                                {
                                    if (kst - diff > Math.Max(MIN_PARAMETER_VALUE, kbt))
                                    {
                                        UpdateKbi(kst - diff, kbt, koa, trainingPeptides, shift);
                                        tempError = trainingPeptides.Sum(x => x.Error);
                                    }
                                }

                                if (tempError > updatedError)
                                {
                                    diff = shift;
                                }
                                else
                                {
                                    updatedError = tempError;
                                }
                            }
                            else
                            {
                                diff = shift;
                            }

                            if (increaseParameter)
                            {
                                variables[0] += diff;
                            }
                            else
                            {
                                variables[0] -= diff;
                            }

                            if (updatedError < bestError)
                            {
                                bestError = updatedError;
                                bestVariables[0] = variables[0];
                                bestVariables[1] = kbt;
                                bestVariables[2] = koa;
                            }
                        }

                        //KBT
                        if (kbt + shift < Math.Min(MAX_KBT_VALUE, kst))
                        {
                            UpdateKbi(kst, kbt + shift, koa, trainingPeptides, shift);
                            updatedError = trainingPeptides.Sum(x => x.Error);
                            increaseParameter = true;
                        }
                        else
                        {
                            updatedError = double.PositiveInfinity;
                        }
                        if (!(updatedError < previousError) && kbt - shift > MIN_PARAMETER_VALUE)
                        {
                            UpdateKbi(kst, kbt - shift, koa, trainingPeptides, shift);
                            updatedError = trainingPeptides.Sum(x => x.Error);
                            increaseParameter = false;
                        }

                        if (previousError > updatedError)
                        {
                            optimizing = true;
                            double diff = (previousError - updatedError) / shift;
                            if (diff > shift * 3)
                            {
                                diff = Math.Round(diff / shift) * shift;
                                double tempError = double.PositiveInfinity;
                                if (increaseParameter)
                                {
                                    if (kbt + diff < Math.Min(MAX_KBT_VALUE, kst))
                                    {
                                        UpdateKbi(kst, kbt + diff, koa, trainingPeptides, shift);
                                        tempError = trainingPeptides.Sum(x => x.Error);
                                    }
                                }
                                else
                                {
                                    if (kbt - diff > MIN_PARAMETER_VALUE)
                                    {
                                        UpdateKbi(kst, kbt - diff, koa, trainingPeptides, shift);
                                        tempError = trainingPeptides.Sum(x => x.Error);
                                    }
                                }

                                if (tempError > updatedError)
                                {
                                    diff = shift;
                                }
                                else
                                {
                                    updatedError = tempError;
                                }
                            }
                            else
                            {
                                diff = shift;
                            }

                            if (increaseParameter)
                            {
                                variables[1] += diff;
                            }
                            else
                            {
                                variables[1] -= diff;
                            }

                            if (updatedError < bestError)
                            {
                                bestError = updatedError;
                                bestVariables[0] = kst;
                                bestVariables[1] = variables[1];
                                bestVariables[2] = koa;
                            }
                        }

                        //KOA                    }
                        if (koa + shift < MAX_KOA_VALUE)
                        {
                            UpdateKbi(kst, kbt, koa + shift, trainingPeptides, shift);
                            updatedError = trainingPeptides.Sum(x => x.Error);
                            increaseParameter = true;
                        }
                        else
                        {
                            updatedError = double.PositiveInfinity;
                        }
                        if (!(updatedError < previousError) && koa - shift > Math.Max(MIN_PARAMETER_VALUE, kst))
                        {
                            UpdateKbi(kst, kbt, koa - shift, trainingPeptides, shift);
                            updatedError = trainingPeptides.Sum(x => x.Error);
                            increaseParameter = false;
                        }

                        if (previousError > updatedError)
                        {
                            optimizing = true;
                            double diff = (previousError - updatedError) / shift;
                            if (diff > shift * 3)
                            {
                                diff = Math.Round(diff / shift) * shift;
                                double tempError = double.PositiveInfinity;
                                if (increaseParameter)
                                {
                                    if (koa + diff < MAX_KOA_VALUE)
                                    {
                                        UpdateKbi(kst, kbt, koa + diff, trainingPeptides, shift);
                                        tempError = trainingPeptides.Sum(x => x.Error);
                                    }
                                }
                                else
                                {
                                    if (koa - diff > Math.Max(MIN_PARAMETER_VALUE, kst))
                                    {
                                        UpdateKbi(kst, kbt, koa - diff, trainingPeptides, shift);
                                        tempError = trainingPeptides.Sum(x => x.Error);
                                    }
                                }

                                if (tempError > updatedError)
                                {
                                    diff = shift;
                                }
                                else
                                {
                                    updatedError = tempError;
                                }
                            }
                            else
                            {
                                diff = shift;
                            }

                            if (increaseParameter)
                            {
                                variables[2] += diff;
                            }
                            else
                            {
                                variables[2] -= diff;
                            }

                            if (updatedError < bestError)
                            {
                                bestError = updatedError;
                                bestVariables[0] = kst;
                                bestVariables[1] = kbt;
                                bestVariables[2] = variables[2];
                            }
                        }

                        //test the new gradient to ensure we don't get trapped
                        UpdateKbi(variables[0], variables[1], variables[2], trainingPeptides, shift);
                        updatedError = trainingPeptides.Sum(x => x.Error);
                        if (updatedError < bestError)
                        {
                            bestError = updatedError;
                            for (int v = 0; v < variables.Length; v++)
                            {
                                bestVariables[v] = variables[v];
                            }
                        }
                        else //reset to the old
                        {
                            UpdateKbi(bestVariables[0], bestVariables[1], bestVariables[2], trainingPeptides, shift);
                        }

                        previousError = bestError;
                    }

                    //save the fits
                    bestKst = bestVariables[0];
                    bestKbt = bestVariables[1];
                    bestKoa = bestVariables[2];
                    kstList[bestIndex] = bestKst;
                    kbtList[bestIndex] = bestKbt;
                    koaList[bestIndex] = bestKoa;
                }
                List<string> paramsForOutput = new List<string> { "kst\tkbt\tkoa", (bestKst.ToString() + '\t' + bestKbt.ToString() + '\t' + bestKoa.ToString()) };
                File.WriteAllLines(pathForFittedParameters, paramsForOutput);
            }

            //Foreach peptide, apply the kst, kbt, and koa, but optimize for the kbi
            UpdateKbi(bestKst, bestKbt, bestKoa, peptides, ITERATIVE_SHIFT);
            //optimize for best fit
            UpdateKbi(bestKst, bestKbt, bestKoa, peptides, ITERATIVE_SHIFT / 10);
            UpdateKbi(bestKst, bestKbt, bestKoa, peptides, ITERATIVE_SHIFT / 100);

            //use the monte carlo method to estimate the 95% confidence interval
            Parallel.ForEach(Partitioner.Create(0, peptides.Count), new ParallelOptions { MaxDegreeOfParallelism = 8 },
    (partitionRange, loopState) =>
    {
        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
        {
            UpdateKbiConfidenceInterval(bestKst, bestKbt, bestKoa, peptides[i], ITERATIVE_SHIFT);
        }
    });

            List<string> linesToWrite = new List<string>();
            linesToWrite.Add("Peptide\tProtein\tHalf Life\tLowerConfidenceInterval\tUpperConfidenceInterval\tError (MSE)\tSummed Intensity\tNumber of Ratios");
            foreach (PeptideTurnoverObject peptide in peptides)
            {
                linesToWrite.Add(peptide.FullSequence + '\t' + peptide.Protein + '\t' + (Math.Log(2,Math.E)/peptide.Kbi).ToString() + '\t' + (Math.Log(2, Math.E) / peptide.HighKbi).ToString() + '\t' + (Math.Log(2, Math.E) / peptide.LowKbi).ToString() + '\t' + peptide.Error.ToString() + '\t' + peptide.TotalIntensity.ToString() + '\t' + peptide.Timepoints.Length.ToString());
            }
            //output each peptide with its sequence, kbi, 95% confidence interval, and protein
            File.WriteAllLines(Path.Combine(directory, filename + "_TurnoverResults.txt"), linesToWrite);
            return new PoolParameters(bestKst, bestKbt, bestKoa);
        }
        

        public static void OptimizeFit(double[] bestVariables, double kst, double kbt, double koa, PeptideTurnoverObject peptide )
        {
            bool optimizing = true;
            double updatedError = 0;
            bool increaseParameter = true;
            UpdateKbi(kst, kbt, koa, peptide, ITERATIVE_SHIFT);
            double previousError = peptide.Error;
            double bestError = previousError;

            while (optimizing)
            {
                optimizing = false;
                kst = bestVariables[0];
                kbt = bestVariables[1];
                koa = bestVariables[2];
                double[] variables = new double[] { kst, kbt, koa };
                
                //KST
                if (kst + ITERATIVE_SHIFT < Math.Min(MAX_KST_VALUE, koa))
                {
                    UpdateKbi(kst + ITERATIVE_SHIFT, kbt, koa, peptide, ITERATIVE_SHIFT);
                    updatedError = peptide.Error;
                    increaseParameter = true;
                }
                else
                {
                    updatedError = double.PositiveInfinity;
                }
                if (!(updatedError < previousError) && kst - ITERATIVE_SHIFT > Math.Max(MIN_PARAMETER_VALUE, kbt))
                {
                    UpdateKbi(kst - ITERATIVE_SHIFT, kbt, koa, peptide, ITERATIVE_SHIFT);
                    updatedError = peptide.Error;
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
                            if (kst + diff < Math.Min(MAX_KST_VALUE, koa))
                            {
                                UpdateKbi(kst + diff, kbt, koa, peptide, ITERATIVE_SHIFT);
                                tempError = peptide.Error;
                            }
                        }
                        else
                        {
                            if (kst - diff > Math.Max(MIN_PARAMETER_VALUE, kbt))
                            {
                                UpdateKbi(kst - diff, kbt, koa, peptide, ITERATIVE_SHIFT);
                                tempError = peptide.Error;
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
                        variables[0] += diff;
                    }
                    else
                    {
                        variables[0] -= diff;
                    }

                    if (updatedError < bestError)
                    {
                        bestError = updatedError;
                        bestVariables[0] = variables[0];
                        bestVariables[1] = kbt;
                        bestVariables[2] = koa;
                    }
                }

                //KBT
                if (kbt + ITERATIVE_SHIFT < Math.Min(MAX_KBT_VALUE, kst))
                {
                    UpdateKbi(kst, kbt + ITERATIVE_SHIFT, koa, peptide, ITERATIVE_SHIFT);
                    updatedError = peptide.Error;
                    increaseParameter = true;
                }
                else
                {
                    updatedError = double.PositiveInfinity;
                }
                if (!(updatedError < previousError) && kbt - ITERATIVE_SHIFT > MIN_PARAMETER_VALUE)
                {
                    UpdateKbi(kst, kbt - ITERATIVE_SHIFT, koa, peptide, ITERATIVE_SHIFT);
                    updatedError = peptide.Error;
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
                            if (kbt + diff < Math.Min(MAX_KBT_VALUE, kst))
                            {
                                UpdateKbi(kst, kbt + diff, koa, peptide, ITERATIVE_SHIFT);
                                tempError = peptide.Error;
                            }
                        }
                        else
                        {
                            if (kbt - diff > MIN_PARAMETER_VALUE)
                            {
                                UpdateKbi(kst, kbt - diff, koa, peptide, ITERATIVE_SHIFT);
                                tempError = peptide.Error;
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
                        variables[1] += diff;
                    }
                    else
                    {
                        variables[1] -= diff;
                    }

                    if (updatedError < bestError)
                    {
                        bestError = updatedError;
                        bestVariables[0] = kst;
                        bestVariables[1] = variables[1];
                        bestVariables[2] = koa;
                    }
                }

                //KOA                    }
                if (koa + ITERATIVE_SHIFT < MAX_KOA_VALUE)
                {
                    UpdateKbi(kst, kbt, koa + ITERATIVE_SHIFT, peptide, ITERATIVE_SHIFT);
                    updatedError = peptide.Error;
                    increaseParameter = true;
                }
                else
                {
                    updatedError = double.PositiveInfinity;
                }
                if (!(updatedError < previousError) && koa - ITERATIVE_SHIFT > Math.Max(MIN_PARAMETER_VALUE, kst))
                {
                    UpdateKbi(kst, kbt, koa - ITERATIVE_SHIFT, peptide, ITERATIVE_SHIFT);
                    updatedError = peptide.Error;
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
                            if (koa + diff < MAX_KOA_VALUE)
                            {
                                UpdateKbi(kst, kbt, koa + diff, peptide, ITERATIVE_SHIFT);
                                tempError = peptide.Error;
                            }
                        }
                        else
                        {
                            if (koa - diff > Math.Max(MIN_PARAMETER_VALUE, kst))
                            {
                                UpdateKbi(kst, kbt, koa - diff, peptide, ITERATIVE_SHIFT);
                                tempError = peptide.Error;
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
                        variables[2] += diff;
                    }
                    else
                    {
                        variables[2] -= diff;
                    }

                    if (updatedError < bestError)
                    {
                        bestError = updatedError;
                        bestVariables[0] = kst;
                        bestVariables[1] = kbt;
                        bestVariables[2] = variables[2];
                    }
                }

                //test the new gradient to ensure we don't get trapped
                UpdateKbi(variables[0], variables[1], variables[2], peptide, ITERATIVE_SHIFT);
                updatedError = peptide.Error;
                if (updatedError < bestError)
                {
                    bestError = updatedError;
                    for (int v = 0; v < variables.Length; v++)
                    {
                        bestVariables[v] = variables[v];
                    }
                }
                else //reset to the old
                {
                    UpdateKbi(bestVariables[0], bestVariables[1], bestVariables[2], peptide, ITERATIVE_SHIFT);
                }

                previousError = bestError;
            }
        }

        public static void UpdateKbi(double kst, double kbt, double koa, List<PeptideTurnoverObject> peptides, double ITERATIVE_SHIFT)
        {
            Parallel.ForEach(Partitioner.Create(0, peptides.Count), new ParallelOptions { MaxDegreeOfParallelism = 8 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        UpdateKbi(kst, kbt, koa, peptides[i], ITERATIVE_SHIFT);
                    }
                });
        }

        private static void UpdateKbi(double kst, double kbt, double koa, PeptideTurnoverObject peptide, double ITERATIVE_SHIFT)
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
                peptide.Error = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi, timepoints, relativeFractions);
                peptide.TemporaryError = peptide.Error;
            }
            double originalError = peptide.Error;

            double updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi + ITERATIVE_SHIFT, timepoints, relativeFractions);
            bool increaseKbi = true;
            if (!(updatedError < originalError))
            {
                updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi - ITERATIVE_SHIFT, timepoints, relativeFractions);
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
                        tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi + diff, timepoints, relativeFractions);
                    }
                    else
                    {
                        if (kbi - diff < MIN_KBI)
                        {
                            diff = kbi - MIN_KBI - ITERATIVE_SHIFT;
                        }
                        tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi - diff, timepoints, relativeFractions);
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

                    updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi + ITERATIVE_SHIFT, timepoints, relativeFractions);
                    increaseKbi = true;
                    if (!(updatedError < originalError))
                    {
                        updatedError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi - ITERATIVE_SHIFT, timepoints, relativeFractions);
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
                            tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi + diff, timepoints, relativeFractions);
                        }
                        else
                        {
                            if (kbi - diff < MIN_KBI)
                            {
                                diff = kbi - MIN_KBI - ITERATIVE_SHIFT;
                            }
                            tempError = CalculateErrorForThreeCompartmentModelFit(kst, kbt, koa, kbi - diff, timepoints, relativeFractions);
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
        /// <param name="koa">Constant for amino acid compartment to waste</param>
        /// <param name="kbi">Constant for the protein of interest's degredation</param>
        private static double CalculateErrorForThreeCompartmentModelFit(double kst, double kbt, double koa, double kbi, double[] timepoints, double[] relativeFraction)
        {
            double[] residuals = new double[timepoints.Length];
            double[] predictedValues = PredictRelativeFractionUsingThreeCompartmentModel(kst, kbt, koa, kbi, timepoints);

            for (int i = 0; i < timepoints.Length; i++)
            {
                residuals[i] = relativeFraction[i] - predictedValues[i];
            }

            //SSE is the sum, this is the average
            return residuals.Average(x => Math.Pow(x, 2));
        }

        //assumes timepoints are sorted
        public static double[] PredictRelativeFractionUsingThreeCompartmentModel(double kst, double kbt, double koa, double kbi, double[] timepoints)
        {
            double sumOfStOaBt = kst + kbt + koa;
            double sqrtForUV = Math.Sqrt(Math.Pow(sumOfStOaBt, 2) - 4 * koa * kbt);
            double u = (sumOfStOaBt - sqrtForUV) / 2;
            double v = (sumOfStOaBt + sqrtForUV) / 2;
            double yu = koa * kbi * (u - kbt) / ((u - v) * (u - kbi) * u);
            double yv = koa * kbi * (v - kbt) / ((v - u) * (v - kbi) * v);
            double ykbi = koa * (kbi - kbt) / ((u - kbi) * (v - kbi));

            double[] predictedValues = new double[timepoints.Length];
            double previousTimepoint = double.NegativeInfinity;
            for (int i = 0; i < timepoints.Length; i++)
            {
                double t = timepoints[i];
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


        public static void UpdateKbiConfidenceInterval(double kst, double kbt, double koa, PeptideTurnoverObject peptide, double ITERATIVE_SHIFT)
        {
            double[] timepoints = peptide.Timepoints;
            //int numMeasurements = timepoints.Length;
            double[] relativeFraction = peptide.RelativeFractions;
            //double[] residuals = new double[numMeasurements];
            double originalKbi = peptide.Kbi; //save it so it doesn't get overwritten
            //double[] distinctPredictedValues = PredictRelativeFractionUsingThreeCompartmentModel(kst, kbt, koa, peptide.Kbi, timepoints).Distinct().ToArray();
            //double[] uniqueTimePoints = timepoints.Distinct().ToArray();
            List<List<double>> relativeFractionsForEachTimePoint = new List<List<double>>();
            double previousTimepoint = timepoints[0];

            int previousI = 0;
            //get predicted values
            double[] predicted = PredictRelativeFractionUsingThreeCompartmentModel(kst, kbt, koa, peptide.Kbi, timepoints);
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
                        sigmasForEachTimePoint.Add(Math.Sqrt(sigmasForThisTimepoint.Sum() / (Math.Max(sigmasForThisTimepoint.Count - 1,1)))); //prevent infinity if only one point
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
            int lowPercentile = (int)Math.Round(NUM_SIMULATIONS * 2.5 / 100); //2.5 percentile
            int highPercentile = (int)Math.Round(NUM_SIMULATIONS * 97.5 / 100); //97.5 percentile

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
                    if(simulatedValue>1)
                    {
                        simulatedValue = 1;
                    }
                    if(simulatedValue<0)
                    {
                        simulatedValue = 0;
                    }
                    simulatedData[i] = simulatedValue;
                    //simulatedData[i] = distinctPredictedValues[i] + NormInv(rng.NextDouble(), 0, sigma);
                }
                peptide.RelativeFractions = simulatedData;

                //optimize on the simulated data
                UpdateKbi(kst, kbt, koa, peptide, ITERATIVE_SHIFT);
                UpdateKbi(kst, kbt, koa, peptide, ITERATIVE_SHIFT/10);
                UpdateKbi(kst, kbt, koa, peptide, ITERATIVE_SHIFT/100);
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
            //    UpdateKbi(kst, kbt, koa, peptide, ITERATIVE_SHIFT);
            //    bootstrapKbis[s] = peptide.Kbi;
            //}
            peptide.Kbi = originalKbi;
            peptide.RelativeFractions = relativeFraction;
            peptide.Timepoints = timepoints;
            peptide.Error = originalError;
            Array.Sort(bootstrapKbis);
            peptide.LowKbi = bootstrapKbis[lowPercentile]; //2.5 percentile
            peptide.HighKbi = bootstrapKbis[highPercentile]; //97.5 percentile
        }


        /// <summary>
        /// Given a probability, a mean, and a standard deviation, an x value can be calculated.
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
            // using Epsilon is wrong; see link above for reference to 0.02425 value
            //const double pLow = double.Epsilon;
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
