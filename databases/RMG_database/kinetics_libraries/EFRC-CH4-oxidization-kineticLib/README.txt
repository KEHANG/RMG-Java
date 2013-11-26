This is a Kinetic Lib from ERFC ch4 oxidization, but takes out the activated radicals like CH* and OH* and related reactions.
Eliminate AR, HE, N2 species.

After using this seed mechanism and generating the RMG mechanism, you should do the following thins:
1. open chem.inp;
2. knock out the Ne in ELEMENT, SPECIES, THERMO;
3. add CH* and OH* to SPECIES, THERMO, add CH*, OH*, Ar, He and N2 to REACTIONS from retrieveInfo/retrieveRxns.txt and retrieveInfo/retrieveThermo.txt;
4. open tran.dat;
5. add CH* and OH* to it from retrieveInfo/retrieveTrans.txt;
