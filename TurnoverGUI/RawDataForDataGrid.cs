using System.IO;

namespace AppleTurnover
{
    public class RawDataForDataGrid
    {
        public string FilePath { get; private set; }
        public string Filename { get; private set; }
        public RawDataForDataGrid(string path)
        {
            FilePath = path;
            Filename = Path.GetFileName(path);
        }
    }
}
