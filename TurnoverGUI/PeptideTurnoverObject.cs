using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TurnoverGUI
{
    public class PeptideTurnoverObject
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public string Protein { get; set; }
        public double[] Timepoints { get; set; }
        public double[] RelativeFractions { get; set; }
        public string[] Filenames { get; set; }
        public double[] Intensities { get; set; }
        public double Kbi { get; set; }
        public double Error { get; set; }
        public double TemporaryError { get; set; }
        public double TotalIntensity { get; set; }
        public double LowKbi { get; set; }
        public double HighKbi { get; set; }
        public int StartResidue { get; set; }
        public int EndResidue { get; set; }
        public Dictionary<int, string> ModDictionary { get; set; }
        public string ErrorString { get; private set; }
        public string FileName { get; set; }

        public PeptideTurnoverObject(string fullSequence, double[] timepoints, double[] rFValues, string[] filenames, double[] intensities, double totalIntensity, string fileName, string protein = "")
        {
            FullSequence = fullSequence;
            BaseSequence = CleanSeq(fullSequence);
            Protein = protein;

            //sort by timepoints
            for (int i = 1; i < timepoints.Length; i++)
            {
                int j = i;
                while (j > 0 && timepoints[j] < timepoints[j - 1])
                {
                    double temp = timepoints[j - 1];
                    timepoints[j - 1] = timepoints[j];
                    timepoints[j] = temp;

                    temp = rFValues[j - 1];
                    rFValues[j - 1] = rFValues[j];
                    rFValues[j] = temp;

                    string tempString = filenames[j - 1];
                    filenames[j - 1] = filenames[j];
                    filenames[j] = tempString;

                    double tempIntensity = intensities[j - 1];
                    intensities[j - 1] = intensities[j];
                    intensities[j] = tempIntensity;
                    j--;
                }
            }

            //Array.Sort(timepoints, rFValues);
            Timepoints = timepoints;
            RelativeFractions = rFValues;
            Filenames = filenames;
            Intensities = intensities;
            Error = double.PositiveInfinity;
            TemporaryError = Error;
            TotalIntensity = totalIntensity;
            Kbi = 0.4; //arbitrary starting value
            LowKbi = 0;
            HighKbi = 2;
            ModDictionary = new Dictionary<int, string>();
            FileName = fileName;
        }

        public PeptideTurnoverObject Copy(string protein)
        {
            return new PeptideTurnoverObject(FullSequence, Timepoints, RelativeFractions, Filenames, Intensities, TotalIntensity, FileName, protein);
        }

        public void UpdateError()
        {
            Error = TemporaryError;
        }

        //removes all modifications (in parenthesis) and returns the base sequence
        public static string CleanSeq(string seq)
        {
            bool ModificationOn = false;
            string ModificationName = "";
            string cleanedSeq = "";
            char[] seqArray = seq.ToCharArray();
            foreach (char amino_acid in seqArray) //if there are synonymous peaks, then the sequences must be identical or possess ambiguities that will be caught later
            {
                if (amino_acid == ']') //only occurs at end of mod
                {
                    ModificationOn = false;
                }
                if (ModificationOn == true) //only occurs if "(" already found
                {
                    ModificationName += amino_acid;
                }
                if (amino_acid == '[') //start collecting PTM name
                {
                    ModificationOn = true;
                }
                if (ModificationOn == false && amino_acid != ']')
                {
                    cleanedSeq += amino_acid;
                }
            }
            return cleanedSeq;
        }

        public void SetErrorString()
        {
            double test = Math.Round(Error, 6);
            ErrorString = test.ToString();
        }
    }
}
