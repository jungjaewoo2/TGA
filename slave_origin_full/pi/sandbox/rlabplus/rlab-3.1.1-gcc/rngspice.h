// Copyright (C) 2003-2013 Marijan Kostrun
//   part of rlabplus for linux project on rlabplus.sourceforge.net
//
// ngspice for rlabplus
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ./COPYING

extern Ent *  ent_ngspice_Initialize(int nargs, Datum args[]);
extern Ent *  ent_ngspice_Finalize (int nargs, Datum args[]);
extern Ent *  ent_ngspice_Command (int nargs, Datum args[]);
extern Ent *  ent_ngspice_RunCkt (int nargs, Datum args[]);
extern Ent *  ent_ngspice_Circuit (int nargs, Datum args[]);
extern Ent *  ent_ngspice_GetVals (int nargs, Datum args[]);
extern Ent *  ent_ngspice_IsRunning (int nargs, Datum args[]);

Bltin rlab_ngspice_bltin[] = {
  {BLTIN, "init",      ent_ngspice_Initialize},
  {BLTIN, "kill",      ent_ngspice_Finalize},
  {BLTIN, "cmd",       ent_ngspice_Command},
  {BLTIN, "runckt",    ent_ngspice_RunCkt},
  {BLTIN, "circuit",   ent_ngspice_Circuit},
  {BLTIN, "getvals",   ent_ngspice_GetVals},
  {BLTIN, "isrunning", ent_ngspice_IsRunning},
  {0, 0, 0},
};
