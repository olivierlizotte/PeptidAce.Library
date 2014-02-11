using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce.Utilities.Interfaces
{
    public delegate void UpdateMethod(string message);

    public class ConSolBasic : IConSol
    {
        protected UpdateMethod funcUpdate = null;

        /// <summary>
        /// Store the list of online users. Wish I had a ConcurrentList. 
        /// </summary>
        protected Logger Logs;
        public ConSolBasic(UpdateMethod updateMethod = null)
        {
            Logs = new Logger();
            funcUpdate = updateMethod;
        }

        public virtual void WriteLine(string msg)
        {
            if (funcUpdate != null)
                funcUpdate(msg);
            else
                Console.WriteLine(msg);
        }

        public void UpdateLogFile()
        {
            Logs.UpdateLogFile();
        }
    }

    public class Logger
    {
        public static string FOLDER = "";
        public static object synchroton = new object();

        private List<string> NewLines = new List<string>(100);

        public Logger()
        {
            System.Threading.TimerCallback tc = new System.Threading.TimerCallback(Tick);
            System.Threading.Timer objTimer = new System.Threading.Timer(tc);
            objTimer.Change(0, 150000);
        }

        //Update the log file
        private void Tick(object data)
        {
            UpdateLogFile();
        }

        public void Add(string msg, string username)
        {
            lock (synchroton)
            {
                NewLines.Add(username + " : " + msg);
            }
        }


        public void UpdateLogFile()
        {
            try
            {
                lock (synchroton)
                {
                    if (NewLines.Count > 0)
                    {
                        //lets create a log file with the name of the user. (needs write access to the executing folder)
                        if (string.IsNullOrEmpty(FOLDER))
                        {
                            FOLDER = AppDomain.CurrentDomain.BaseDirectory + "Logs" + System.IO.Path.DirectorySeparatorChar;
                            if (!System.IO.Directory.Exists(FOLDER))
                                System.IO.Directory.CreateDirectory(FOLDER);
                        }

                        string fileName = System.DateTime.Today.Day + "_" + System.DateTime.Today.Month + "_" + System.DateTime.Today.Year + ".txt";
                        System.IO.StreamWriter outFile = new System.IO.StreamWriter(FOLDER + System.IO.Path.DirectorySeparatorChar + fileName, true);

                        foreach (string line in NewLines)
                            outFile.WriteLine(line.ToString());

                        outFile.Close();
                        NewLines.Clear();
                    }
                }
            }
            catch (System.Exception)
            {
                //Exceptionnally here, since its a log failing, do nothing!
            }
        }
    }
}
