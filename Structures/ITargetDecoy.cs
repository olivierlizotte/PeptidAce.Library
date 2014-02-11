using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PeptidAce
{
    public interface ITargetDecoy
    {
        bool Target { get; }
        bool Decoy { get; }
    }
}
