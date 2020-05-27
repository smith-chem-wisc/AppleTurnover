using System;
using System.Collections.Generic;
using System.Text;

namespace AppleTurnover
{
    public class Settings
    {
        public int MinValidValuesTotal { get; set; }
        public int MinValidValuesPerTimepoint { get; set; }
        public bool UseBadRatios { get; set; }
        public bool RemoveMessyPeptides { get; set; }
        public SearchEngine UpstreamProgram {get;set;}
        public Settings(int minValidValuesTotal = 6, int minValidValuesPerTimepoint= 3, bool useBadRatios = false, SearchEngine engine = SearchEngine.MetaMorpheus, bool removeMessyPeptides = true)
        {
            MinValidValuesTotal = minValidValuesTotal;
            MinValidValuesPerTimepoint = minValidValuesPerTimepoint;
            UseBadRatios = useBadRatios;
            UpstreamProgram = engine;
            RemoveMessyPeptides = removeMessyPeptides;
        }

        public enum SearchEngine
        {
            MaxQuant,
            MetaMorpheus
        }
    }
}
