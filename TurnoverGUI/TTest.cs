using FlashLFQ;
using Meta.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.ComTypes;
using System.Text;

namespace AppleTurnover
{
    public static class TTest
    {
        public static void CompareProteinsAcrossFiles(List<string> filenames, List<PeptideTurnoverObject> allProteins, Dictionary<string, PoolParameters> poolParameterDictionary)
        {
            if (filenames.Count < 2)
            {
                return;
            }
            string outputDirectory = Path.GetDirectoryName(filenames.First());
            Directory.CreateDirectory(Path.Combine(outputDirectory, "StatisticalComparisons"));
            for (int i = 0; i < filenames.Count; i++)
            {
                string fileOne = filenames[i];
                PoolParameters paramsOne = poolParameterDictionary[fileOne];
                //get the proteins for this file
                List<PeptideTurnoverObject> proteinsForFileOne = allProteins.Where(x => fileOne.Equals(x.FileName)).OrderBy(x => x.Protein).ToList();

                for (int j = i + 1; j < filenames.Count; j++)
                {
                    string fileTwo = filenames[j];
                    PoolParameters paramsTwo = poolParameterDictionary[fileTwo];
                    List<string> linesToWrite = new List<string>();
                    //add header
                    linesToWrite.Add("Protein\tFold Change\tNeg. log(p-Value)\tHalf-life " + fileOne + "\tHalf-life " + fileTwo);

                    List<PeptideTurnoverObject> proteinsForFileTwo = allProteins.Where(x => fileTwo.Equals(x.FileName)).OrderBy(x => x.Protein).ToList();

                    //get the overlap between them
                    int a = 0;
                    int b = 0;
                    while (a < proteinsForFileOne.Count && b < proteinsForFileTwo.Count)
                    {
                        PeptideTurnoverObject proteinOne = proteinsForFileOne[a];
                        PeptideTurnoverObject proteinTwo = proteinsForFileTwo[b];
                        int comparison = (proteinOne.Protein).CompareTo(proteinTwo.Protein);
                        if (comparison == 0)
                        {
                            //do the comparison (t-test of montecarlos, which dramatically overestimates the sample size)
                            //Sample sampleOne = new Sample(proteinOne.MonteCarloKbis.Select(x => Math.Log10(2) / x));
                            //Sample sampleTwo = new Sample(proteinTwo.MonteCarloKbis.Select(x => Math.Log10(2) / x));
                            //TestResult result = Sample.StudentTTest(sampleOne, sampleTwo);
                            //linesToWrite.Add(proteinOne.Protein + "\t" + (Math.Log2(sampleTwo.Median) - Math.Log2(sampleOne.Median)).ToString() + '\t' +
                            //    (-1*Math.Log10(result.Probability)).ToString() + '\t' + (Math.Log10(2) / proteinOne.Kbi).ToString() + '\t' + (Math.Log10(2) / proteinTwo.Kbi).ToString());

                            //do the comparison (t-test of normalized ratios for all timepoints)
                            double averageKbi = (proteinOne.Kbi + proteinTwo.Kbi) / 2;
                            double normalizedHalfLife = Math.Log(2) / (averageKbi); //this is the day we're going to normalize all of the relative fractions to

                            //create an array of a single value (the normalized timepoint) to create a new timepoint array
                            double[] comparisonTimepointsOne = new double[proteinOne.Timepoints.Length];
                            double[] comparisonTimepointsTwo = new double[proteinTwo.Timepoints.Length];
                            for (int index = 0; index < comparisonTimepointsOne.Length; index++)
                            {
                                comparisonTimepointsOne[index] = normalizedHalfLife;
                            }
                            for (int index = 0; index < comparisonTimepointsTwo.Length; index++)
                            {
                                comparisonTimepointsTwo[index] = normalizedHalfLife;
                            }

                            //predict the expected values for the ratios of protein one based on the fit of the comparison
                            double[] expectedOriginalRatiosOne = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(paramsOne.Kst, paramsOne.Kbt, paramsOne.Kao, averageKbi, proteinOne.Timepoints);
                            //predict the expected values for the ratios of proteoform one based on the fit of the comparison if they were all at the same normalized timepoint
                            double[] expectedUpdatedRatiosOne = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(paramsOne.Kst, paramsOne.Kbt, paramsOne.Kao, averageKbi, comparisonTimepointsOne);
                            
                            //do the same thing with protein two
                            double[] expectedOriginalRatiosTwo = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(paramsTwo.Kst, paramsTwo.Kbt, paramsTwo.Kao, averageKbi, proteinTwo.Timepoints);
                            double[] expectedUpdatedRatiosTwo = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(paramsTwo.Kst, paramsTwo.Kbt, paramsTwo.Kao, averageKbi, comparisonTimepointsTwo);

                            //create empty arrays for the normalized ratios
                            double[] normalizedRatiosOne = new double[expectedOriginalRatiosOne.Length];
                            double[] normalizedRatiosTwo = new double[expectedOriginalRatiosTwo.Length];

                            //calculate the normalized ratios by subtracting the expected ratio (so that we are measuring the residual between the point and the comparison fit) and then adding the ratio of the comparison fit at the normalized timepoint. 
                            for (int index = 0; index < proteinOne.RelativeFractions.Length; index++)
                            {
                                //the normalized ratio is equal to the original ratio minus the original fit to the data plus the fit if the kbi was averaged
                                normalizedRatiosOne[index] = proteinOne.RelativeFractions[index] - expectedOriginalRatiosOne[index] + expectedUpdatedRatiosOne[index];
                            }
                            for (int index = 0; index < proteinTwo.RelativeFractions.Length; index++)
                            {
                                normalizedRatiosTwo[index] = proteinTwo.RelativeFractions[index] - expectedOriginalRatiosTwo[index] + expectedUpdatedRatiosTwo[index];
                            }
                            Sample sampleOne = new Sample(normalizedRatiosOne);
                            Sample sampleTwo = new Sample(normalizedRatiosTwo);
                            TestResult result = Sample.StudentTTest(sampleOne, sampleTwo);
                            linesToWrite.Add(proteinOne.Protein + "\t" + (Math.Log2(Math.Log(2) / proteinTwo.Kbi) - Math.Log2(Math.Log(2) / proteinOne.Kbi)).ToString() + '\t' +
                                (-1 * Math.Log10(result.Probability)).ToString() + '\t' + (Math.Log(2) / proteinOne.Kbi).ToString() + '\t' + (Math.Log(2) / proteinTwo.Kbi).ToString());

                            a++;
                            b++;
                        }
                        else if (comparison < 0)
                        {
                            a++;
                        }
                        else
                        {
                            b++;
                        }
                    }
                    File.WriteAllLines(Path.Combine(outputDirectory, "StatisticalComparisons", "Comparison_" + Path.GetFileNameWithoutExtension(fileOne) + "vs" + Path.GetFileNameWithoutExtension(fileTwo) + ".tsv"), linesToWrite);
                }
            }
        }

        public static void CompareProteoformsWithinFiles(List<string> filenames, List<PeptideTurnoverObject> allProteins, Dictionary<string, PoolParameters> poolParameterDictionary)
        {
            for (int fileIndex = 0; fileIndex < filenames.Count; fileIndex++)
            {
                string filename = filenames[fileIndex];
                PoolParameters poolParams = poolParameterDictionary[filename];
                List<PeptideTurnoverObject> proteinsForThisFile = allProteins.Where(x => filename.Equals(x.FileName)).OrderBy(x => x.Proteoform).ToList();
                List<string> linesToWrite = new List<string>();
                linesToWrite.Add("Proteoform A\tProteoform B\tHalf-life A\tHalf-life B\tLog2(Fold Change)\tNeg. log(p-Value)");

                int indexOfNextProteoformFamily = 0;
                for (int i = 0; i < proteinsForThisFile.Count; i++)
                {
                    string currentProtein = proteinsForThisFile[i].Proteoform.Split('_')[0];

                    //find last index for this proteoform family
                    indexOfNextProteoformFamily++;
                    for (; indexOfNextProteoformFamily < proteinsForThisFile.Count; indexOfNextProteoformFamily++)
                    {
                        if (!currentProtein.Equals(proteinsForThisFile[indexOfNextProteoformFamily].Proteoform.Split('_')[0]))
                        {
                            break;
                        }
                    }

                    for (; i < indexOfNextProteoformFamily; i++)
                    {
                        PeptideTurnoverObject proteinOne = proteinsForThisFile[i];

                        //see if it has a localized mod (or localized unmodified site), otherwise skip
                        string[] proteoformOne = proteinOne.Proteoform.Split('@').ToArray();
                        if (proteoformOne.Length == 2)
                        {
                            for (int j = i + 1; j < indexOfNextProteoformFamily; j++)
                            {
                                PeptideTurnoverObject proteinTwo = proteinsForThisFile[j];
                                string[] proteoformTwo = proteinTwo.Proteoform.Split('@').ToArray();

                                //if these are a pair for the same modification site, then do the comparison
                                if (proteoformTwo.Length == 2 && proteoformOne[1].Equals(proteoformTwo[1]))
                                {
                                    //do the comparison (t-test of normalized ratios for all timepoints)
                                    double averageKbi = (proteinOne.Kbi + proteinTwo.Kbi) / 2;
                                    double normalizedHalfLife = Math.Log(2) / (averageKbi); //this is the day we're going to normalize all of the relative fractions to

                                    //create an array of a single value (the normalized timepoint) to create a new timepoint array
                                    double[] comparisonTimepointsOne = new double[proteinOne.Timepoints.Length];
                                    double[] comparisonTimepointsTwo = new double[proteinTwo.Timepoints.Length];
                                    for (int index = 0; index < comparisonTimepointsOne.Length; index++)
                                    {
                                        comparisonTimepointsOne[index] = normalizedHalfLife;
                                    }
                                    for (int index = 0; index < comparisonTimepointsTwo.Length; index++)
                                    {
                                        comparisonTimepointsTwo[index] = normalizedHalfLife;
                                    }

                                    double[] expectedOriginalRatiosOne = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Kao, averageKbi, proteinOne.Timepoints);
                                    double[] expectedUpdatedRatiosOne = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Kao, averageKbi, comparisonTimepointsOne);
                                    double[] expectedOriginalRatiosTwo = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Kao, averageKbi, proteinTwo.Timepoints);
                                    double[] expectedUpdatedRatiosTwo = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Kao, averageKbi, comparisonTimepointsTwo);
                                    double[] normalizedRatiosOne = new double[expectedOriginalRatiosOne.Length];
                                    double[] normalizedRatiosTwo = new double[expectedOriginalRatiosTwo.Length];
                                    for (int index = 0; index < proteinOne.RelativeFractions.Length; index++)
                                    {
                                        //the normalized ratio is equal to the original ratio minus the original fit to the data plus the fit if the kbi was averaged
                                        normalizedRatiosOne[index] = proteinOne.RelativeFractions[index] - expectedOriginalRatiosOne[index] + expectedUpdatedRatiosOne[index];
                                    }
                                    for (int index = 0; index < proteinTwo.RelativeFractions.Length; index++)
                                    {
                                        normalizedRatiosTwo[index] = proteinTwo.RelativeFractions[index] - expectedOriginalRatiosTwo[index] + expectedUpdatedRatiosTwo[index];
                                    }
                                    Sample sampleOne = new Sample(normalizedRatiosOne);
                                    Sample sampleTwo = new Sample(normalizedRatiosTwo);
                                    TestResult result = Sample.StudentTTest(sampleOne, sampleTwo);

                                    try //sometimes crashes if stdev is zero
                                    {
                                        linesToWrite.Add(proteinOne.Proteoform + "\t" + proteinTwo.Proteoform + '\t' + (Math.Log(2) / proteinOne.Kbi).ToString() + '\t' + (Math.Log(2) / proteinTwo.Kbi).ToString() + '\t' +
                                            (Math.Log2((Math.Log(2) / proteinTwo.Kbi)) - Math.Log2((Math.Log(2) / proteinOne.Kbi))).ToString() + '\t' + (-1 * Math.Log(result.Probability)).ToString());
                                    }
                                    catch
                                    {
                                        linesToWrite.Add(proteinOne.Proteoform + "\t" + proteinTwo.Proteoform + '\t' + (Math.Log(2) / proteinOne.Kbi).ToString() + '\t' + (Math.Log(2) / proteinTwo.Kbi).ToString() + '\t' +
                                        (Math.Log2(sampleTwo.Median) - Math.Log2(sampleOne.Median)).ToString() + '\t' + "NA");
                                    }
                                }
                            }
                        }
                    }
                    i--;
                }

                File.WriteAllLines(Path.Combine(Path.GetDirectoryName(filename), Path.GetFileNameWithoutExtension(filename) + "_Results", Path.GetFileNameWithoutExtension(filename) + "_ProteoformAnalysis.tsv"), linesToWrite);
            }
        }
    }
}