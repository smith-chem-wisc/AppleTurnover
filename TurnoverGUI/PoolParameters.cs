using System;
using System.Collections.Generic;
using System.Text;

namespace AppleTurnover
{
    public class PoolParameters
    {
        public double Kst { get; set; }
        public double Kbt { get; set; }
        public double Koa { get; set; }
        public PoolParameters(double kst, double kbt, double koa)
        {
            Kst = kst;
            Kbt = kbt;
            Koa = koa;
        }
    }
}
