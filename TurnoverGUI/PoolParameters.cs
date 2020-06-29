using System;
using System.Collections.Generic;
using System.Text;

namespace AppleTurnover
{
    public class PoolParameters
    {
        public double Kst { get; set; }
        public double Kbt { get; set; }
        public double Kao { get; set; }
        public double DeltaX { get; set; }
        public PoolParameters(double kst, double kbt, double kao, double deltaX)
        {
            Kst = kst;
            Kbt = kbt;
            Kao = kao;
            DeltaX = deltaX;
        }
    }
}
