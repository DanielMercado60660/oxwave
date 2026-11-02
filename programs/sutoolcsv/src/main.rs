// Copyright (c) 2008, Colorado School of Mines,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or
// without modification, are permitted provided that the following
// conditions are met:
//
//     *  Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     *  Redistributions in binary form must reproduce the above
//        copyright notice, this list of conditions and the following
//        disclaimer in the documentation and/or other materials provided
//        with the distribution.
//     *  Neither the name of the Colorado School of Mines nor the names of
//        its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
//
// Warranty Disclaimer:
// THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Rust port for OxWave, a Rust port of CWP/SU.

//! SUTOOLCSV - Tools for Comma Separated Values and fixed text files.
//! Port of su/main/headers/sutoolcsv.c (Revision 1.2).
//! Credits: Andre Latour. Started from sugeom by Fernando Roxo <petro@roxo.org>.

#![allow(non_snake_case)]
use oxwave_core::cio::fputs_bytes;
use oxwave_core::prelude::*;
mod checks;
mod presets;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                          ",
    " SUTOOLCSV - Tools for Comma Separated Values and fixed text files.       ",
    "                                                                          ",
    " sutoolcsv rfile=inx.txt wfile=outx.txt setid=x                           ",
    "           match=sps2 names=sps2 forms=sps2                               ",
    "                                                                          ",
    " Parameter overview:                                                      ",
    "                                                                          ",
    "       rfile= read text values file                                       ",
    "       wfile= write text values file                                      ",
    "       rtype= type of records in input text file (comma separated, fixed) ",
    "       wtype= type of records to output in wfile (comma separated, fixed) ",
    "       names= assign SU names to text values (with SPS2 and SPS1 options) ",
    "       forms= assign output formats to values (with SPS2 and SPS1 options)",
    "       setid= accept data records based on first characters (X,S,R etc.)  ",
    "       match= SU keys to be used by SUGEOMCSV (traces are not input here).",
    "     process= several options to modify/repair values before output.      ",
    "       width= width of records for fixed output (default is 80).          ",
    "      rdelim= rfile delimiter when rtype is csv (default is comma).       ",
    "  nicerecord= allows skipping bad records at start of some input files    ",
    "  maxrecords= maximum number of records to allocate memory for.           ",
    "    unrepeat= help when duplicate fldr values exist (default is off)      ",
    "      scalco= check size of sx,sy,gx,gy coordinate values                 ",
    "      scalel= check size of elevation and related values                  ",
    "       outid= output data records with these first characters (Q,S, etc.) ",
    "                                                                          ",
    " ***********************************************************              ",
    "   To output this documentation:  sutoolcsv 2> tooldoc.txt                ",
    " ***********************************************************              ",
    "                                                                          ",
    " Typical Usage to check SPS files for problems and output:                ",
    "   sutoolcsv rfile=A.txt wfile=B.txt setid=x                              ",
    "             match=sps2 names=sps2 forms=sps2                             ",
    "                                                                          ",
    " Usage to convert SPS to comma-separated output:                          ",
    "   sutoolcsv rfile=A.txt wfile=B.csv setid=x                              ",
    "             match=sps2 names=sps2 forms=sps2                             ",
    "                                                                          ",
    " Usage to convert that comma-separated SPS back to fixed:                 ",
    "   sutoolcsv rfile=B.csv wfile=C.txt                                      ",
    "                                                                          ",
    " *** An important point here is that the SPS options are written          ",
    " *** to records in the output text file in their expanded versions.       ",
    " *** Their expanded versions have all the complicated specifications      ",
    " *** needed for SPS. If necessary, you can paste those records back       ",
    " *** into your original text file, with whatever modifications you need.  ",
    " *** Then use those modifications in another run of this program by       ",
    " *** NOT specifying setid= or match= names= or forms=.                    ",
    "                                                                          ",
    "                                                                          ",
    " ---------------------------------------------------------------------    ",
    "                                                                          ",
    " Starting with a simple example, consider 3 values in a fixed format file ",
    " (say record id, energy source number, and shot elevation).               ",
    "                                                                          ",
    "       S       11  343                                                    ",
    "       S       42  342                                                    ",
    "       S       25  340                                                    ",
    "       C  some comment                                                    ",
    "       S       45  347                                                    ",
    "                                                                          ",
    "  Specifying this on the command line:                                    ",
    "                                                                          ",
    "  sutoolcsv rfile=in.txt wfile=out.csv setid=S match=es                   ",
    "            names=C_su_id,2_es_10,11_selev_15 forms=c_su_id,%.0f          ",
    "                                                                          ",
    "  Results in this set of data records in the output file:                 ",
    "                                                                          ",
    "       S,11,343                                                           ",
    "       S,42,342                                                           ",
    "       S,25,340                                                           ",
    "       S,45,347                                                           ",
    "                                                                          ",
    "  Note the names= specifications. The 2_es_10 specification means read    ",
    "  the value from fixed locations 2 to 10 (card columns) and treat it like ",
    "  SU key es (energy source number). And 11_selev_15 means read fixed      ",
    "  locations 11 to 15 and treat it like SU key selev (shot elevation).     ",
    "                                                                          ",
    "  The forms= specification indicates how you want the numbers formatted   ",
    "  in the output text file. The setid=s specification says to only read    ",
    "  data from records that start with S. The match=es specification says    ",
    "  how these values are going to be merged when you use SUGEOMCSV (it also ",
    "  tells this program what to error-check and warn-check).                 ",
    "                                                                          ",
    "  You will notice in the out.csv file that the data records are preceded  ",
    "  by C_SU records that contain copies of the options you specified on the ",
    "  command line. If you want to re-input the out.csv file in this program  ",
    "  or in SUGEOMCSV, you do not need to re-specify the command line options.",
    "  The options will be read from C_SU records in the out.csv file. This    ",
    "  functionality is not needed for simple files like this 3 value example. ",
    "  But it becomes important for SPS and other complicated files.           ",
    "                                                                          ",
    "       C_SU_MATCH,es                                                      ",
    "       C_SU_SETID,S                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_ID,%.0f                                                       ",
    "       C_SU_NAMES                                                         ",
    "       C_su_id,2_es_10,11_selev_15                                        ",
    "       S,11,343                                                           ",
    "       S,42,342                                                           ",
    "       S,25,340                                                           ",
    "       S,45,347                                                           ",
    "                                                                          ",
    "  Note that commas are still specified in the C_SU_ parameter records     ",
    "  even when the rest of the rfile is fixed format.                        ",
    "                                                                          ",
    "  --------------                                                          ",
    "                                                                          ",
    " There are 2 Parameters which always require command line specification:  ",
    "                                                                          ",
    "       rfile=  read text values file (fixed or comma separated values)    ",
    "                                                                          ",
    "       wfile=  write text values file (fixed or comma separated values)   ",
    "                                                                          ",
    " All other parameters are either not required or can be contained within  ",
    " the input rfile itself.                                                  ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "       rtype= type of records in input text file. Default is csv if the   ",
    "              file name ends in csv, otherwise defaults to fixed.         ",
    "            =csv     comma separated values                               ",
    "            =fixed   This option is usually required for SPS files        ",
    "                     along with other specifications in names= list       ",
    "                     (see extensive examples below).                      ",
    "                                                                          ",
    "       wtype= type of records to output in wfile. Default is csv if the   ",
    "              file name ends in csv, otherwise defaults to fixed.         ",
    "            =csv     comma separated values                               ",
    "            =fixed    This option is required to output SPS files along   ",
    "                      with leading and trailing integers in names= list   ",
    "                      (see extensive examples below).                     ",
    "            =csvchop   Output the raw fields (rather than their values    ",
    "                       converted to numbers and back to characters).      ",
    "                       This can still be used as input to SpreadSheets    ",
    "                       (if extensive repairs or alterations are needed).  ",
    "                   *** Despite outputting the raw fields, all converting, ",
    "                       storing, sorting, and checking continue as normal  ",
    "                       unless you set maxrecords=-1 or -2.                ",
    "                                                                          ",
    " ----------------------------------------------------------------------   ",
    "                                                                          ",
    " The following 4 parameters must be found on the command line or in       ",
    " their corresponding C_SU records in the rfile.                           ",
    "                                                                          ",
    "       match= any number of SU keys needed to find the exact record in    ",
    "             the rfile text file when SUGEOMCSV is used to update the     ",
    "             actual seismic datafile. No traces are input by SUTOOLCSV    ",
    "             but these keys are used herein for various checking purposes.",
    "             These SU keys must also be in the names= list.               ",
    "             Example on command line:                                     ",
    "             match=fldr,tracf                                             ",
    "       match=SPS2 This is just a standard way to specify the match= list  ",
    "                  for SPS Revison 2 files (see the examples below).       ",
    "                  The setid= option must also be X,S, or R.               ",
    "       match=SPS1 This is just a standard way to specify the match= list  ",
    "                  for SPS Revison 1 files (see the examples below).       ",
    "                  The setid= option must also be X,S, or R.               ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "      setid= is used to accept data records based on their first field.   ",
    "        setid=S     means accept data records if their first field is S   ",
    "                          (any characters allowed, such as R,X,cdp,FRED)  ",
    "                    Note: this value is automatically upper-cased unless  ",
    "                          you surround it by double-quotes.               ",
    "                          So s becomes S unless you use double-quotes.    ",
    "        setid=ANY   means read all records (except those starting C_SU)   ",
    "                    and those records have an id field at front.          ",
    "        setid=NONE  means read all records (except those starting C_SU)   ",
    "                    but those records do not have an id field at front.   ",
    "                    (For csv files this means the field before the first  ",
    "                     comma is a value, not an identifier).                ",
    "             Example on command line:                                     ",
    "             setid=S                                                      ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "       names= is used to assign names to values in rfile text file (you   ",
    "              are telling this program the names of values in text file). ",
    "             For files with comma-separated values, a name must be listed ",
    "             sequentially for each field in the rfile text file.          ",
    "             The names must also include the match= SU keys above.        ",
    "             Note C_su_id means this is field used for record acceptance. ",
    "        ***  Read the note   c_su_id IS SPECIAL   later. ***              ",
    "             Special name: null1 (null and any integer) which means       ",
    "                           do not read/output this field. (You can also   ",
    "                           just put nothing between sequential commas).   ",
    "             Special name: numb1 (numb and any integer) which means       ",
    "                           read/output this field even though it is not   ",
    "                           a SU key.                                      ",
    "               Example on command line:                                   ",
    "               names=C_su_id,cdp,null3,cx,cy,,ce                          ",
    "       names=SPS2 This is just a standard way to specify the names= list  ",
    "                  for SPS Revison 2 files (see the examples below).       ",
    "                  The setid= option must also be X,S, or R.               ",
    "       names=SPS1 This is just a standard way to specify the names= list  ",
    "                  for SPS Revison 1 files (see the examples below).       ",
    "                  The setid= option must also be X,S, or R.               ",
    "       names=SPS2ALL  This is a standard way to specify names= to output  ",
    "                  every field from SPS Revison 2 files (see examples).    ",
    "                  The setid= option must also be X,S, or R.               ",
    "       names=SPS1ALL  This is a standard way to specify names= to output  ",
    "                  every field from SPS Revison 1 files (see examples).    ",
    "                  The setid= option must also be X,S, or R.               ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "       forms= is used to assign format specifiers for how to write values ",
    "              to wfile. These are C language formats, but since all values",
    "              are stored internally as double-precision floating point,   ",
    "              only use formats that make sense for that. If there are     ",
    "              fewer specifiers here than names= then the last specifier   ",
    "              is repeated.                                                ",
    "               Example on command line:                                   ",
    "               forms=c_su_id,%.5f                                         ",
    "       forms=SPS2 This is just a standard way to specify the forms= list  ",
    "                  for SPS Revison 2 files (see the examples below).       ",
    "                  The setid= option must also be X,S, or R.               ",
    "       forms=SPS1 The formats for SPS1 option are the same as SPS2 option.",
    "   *** Note that SPS2 and SPS1 files have some fields that are officially ",
    "   *** defined to be character. This program ignores that and assumes     ",
    "   *** anything you try to read is a number and so the options here use   ",
    "   *** a numeric format for every field (except the setid field).         ",
    "   *** If there really are characters in these fields, see wtype=csvchop  ",
    "   *** and process= options for what you might do to alter them.          ",
    "                                                                          ",
    "   Note that SPS options are written to output file in their expanded     ",
    "   version. You can alter these expanded versions and paste them back     ",
    "   to the input file. You might input a text file with a single record    ",
    "   (starting with setid) just to output expanded versions of SPS options. ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "       process= several options which allow modification of values.       ",
    "                rtype=fixed is required. These options were chosen for    ",
    "                their maximum likelyhood to allow repairs of SPS files    ",
    "                without using a SpreadSheet program. But some files will  ",
    "                still require outputting csv and using a SpreadSheet.     ",
    "                These options all have leading and trailing integers which",
    "                specify the characters (card columns) they apply to.      ",
    "              =nn_blank_mm   means set everything from nn to mm to blank. ",
    "              =nn_zero_mm    means set everything from nn to mm to zeros  ",
    "                             (fill with zeros, not just 1 zero).          ",
    "              =nn_trim_mm    set everything from nn to the first + - or   ",
    "                             0-9 digit to blank. And set everything to    ",
    "                             blank from mm back to the last 0-9 digit.    ",
    "              =nn_trimz_mm   set everything from nn to the first + - or   ",
    "                             0-9 digit to blank. And set everything to    ",
    "                             blank from mm back to the last 0-9 digit.    ",
    "                             And set everything that is not a 0-9 to 0    ",
    "                             between the new trimmed ends (except for     ",
    "                             a leading + or -).                           ",
    "              =nn_sub_qqqq_mm subtract the number qqqq from the value     ",
    "                              between nn and mm.                          ",
    "              =nn_add_qqqq_mm add the number qqqq to the value            ",
    "                              between nn and mm.                          ",
    "              =nn_div_qqqq_mm divide the number qqqq into the value       ",
    "                              between nn and mm.                          ",
    "              =nn_mul_qqqq_mm multiply the number qqqq with the value     ",
    "                              between nn and mm.                          ",
    "        The nn and mm numbers do not have to correspond to anything in    ",
    "        names= list. The math options expect one number within nn to mm.  ",
    "        Note these options are performed in the order you specify them.   ",
    "               Example: 11_blank_14,21_sub_1000_30,21_div_3.2808_30       ",
    "                                                                          ",
    "     *Note*  The process= options occur before wtype=csvchop so output    ",
    "             file for csvchop will still have the process= modifications. ",
    "                                                                          ",
    "     *Note*  The c language routines that read-in numbers will stop       ",
    "             at the first non-numeric character. So you may get away      ",
    "             without as much trimming as initially appears needed.        ",
    "             But other languages and programs may not be so forgiving.    ",
    "                                                                          ",
    " -----------------                                                        ",
    " -----------------                                                        ",
    "                                                                          ",
    "       If match= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_MATCH and reads keys from it. ",
    "             Example within the text file:                                ",
    "             C_SU_MATCH,fldr,tracf                                        ",
    "                                                                          ",
    "       If setid= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_SETID and reads id from it.   ",
    "             Example within the text file:                                ",
    "             C_SU_SETID,S                                                 ",
    "                                                                          ",
    "       If names= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_NAMES and reads names from    ",
    "       the record after the C_SU_NAMES record.                            ",
    "             Example within the text file:                                ",
    "             C_SU_NAMES                                                   ",
    "             C_su_id,cdp,null,cx,cy,null,ce                               ",
    "                                                                          ",
    "       If forms= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_FORMS and reads formats from  ",
    "       the record after the C_SU_FORMS record.                            ",
    "             Example within the text file:                                ",
    "             C_SU_FORMS                                                   ",
    "             C_su_id,%.0f,%.2f                                            ",
    "                                                                          ",
    " Note these C_SU_ parameter records can be in any order within text file  ",
    " but the C_SU_NAMES and C_SU_FORMS records must be followed by their      ",
    " correct records. Note also that this program and SUGEOMCSV know that     ",
    " records starting with C_SU are not data records, and will not try to     ",
    " read data values from them even when setid=ANY or NONE.                  ",
    "                                                                          ",
    " *** The previous records will be written to the output file          *** ",
    " *** unless you specify width= negative.                              *** ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "       width= if wtype=fixed this specifies the width for output wfile.   ",
    "              Default is 80 or maximum trailing integer in names= list    ",
    "              whichever is greater. Error if less than maximum trailing   ",
    "              (and a positive value must be greater or equal to 30).      ",
    "            = negative means suppress output of C_SU_ records but still   ",
    "              use abs(width) as width for output wfile.                   ",
    "            = -1 means suppress output of C_SU_ records but still default ",
    "                 to 80 or maximum trailing integer in names= list if      ",
    "                 wtype=fixed (-1 also allowed for wtype= not fixed).      ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "      rdelim= If rtype=fixed you cannot specify this parameter.           ",
    "              If rtype=csv the default is comma. You can specify any      ",
    "              single character here either by itself or surrounded by     ",
    "              double-quotes (needed because some characters such as       ",
    "              semi-colon may have trouble getting through command line).  ",
    "      *Note*  The output always uses commas (if wtype is not fixed).      ",
    "     **Note** Specifying a blank here usually will not give good results  ",
    "              unless the input rfile has exactly 1 blank between numbers. ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "     unrepeat= The default is not to enable this option.                  ",
    "               This option is general but most likely usefull for X-files ",
    "               where the field record number (fldr) increases but then    ",
    "               re-starts at a lower number. Such as 1->7800 then 5->4000. ",
    "               Normally, the finding-logic in SUTOOLCSV and SUGEOMCSV     ",
    "               would not be able to distinguish the first fldr 5 from     ",
    "               the second 5. (For that situation if you do not use this   ",
    "               option you will most likely get multiple layout segment    ",
    "               and channel range warnings for fldr 5 because the same     ",
    "               channel ranges exist twice for fldr 5).                    ",
    "      unrepeat=1 Read the text file and generate an integer from 1 and    ",
    "                 increment by 1 every time the first match= reverses.     ",
    "                 Typically, the first match= is fldr for X-files so this  ",
    "                 increments +1 when fldr is increasing and then decreases,",
    "                 and also increments +1 if fldr is decreasing and then    ",
    "                 increases. The comparison is done using order of records ",
    "                 as they exist in the text file (before sorting herein).  ",
    "                 In SUGEOMCSV this option generates another incrementing  ",
    "                 integer the same way except using the order of traces.   ",
    "                 These two integers are used to match which (fldr) value  ",
    "                 in the traces belongs to which (fldr) value from X-file. ",
    "      unrepeat=  any other integer works the same as 1 in this program but",
    "                 can be specified here in order to maintain consistant    ",
    "                 parameterization with SUGEOMCSV.                         ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "     nicerecord= record number to start trying to read data values from   ",
    "                 the rfile text file (default is 1). The beginning records",
    "                 of some text files are odd (comments and information).   ",
    "                 When the setid= option is not able to reject them,       ",
    "                 specify a record number here where setid= will work.     ",
    "                 (This program also always knows that records starting    ",
    "                 with C_SU are not data records, and will not try to      ",
    "                 read data values from them even when setid=ANY). But     ",
    "                 it will read C_SU parameter records even if they are     ",
    "                 previous to this nicerecord number.                      ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "     maxrecords= maximum number of records to allocate memory for.        ",
    "                 If not specified, this program reads through the records ",
    "                 once and allocates for the number found. Then reads them ",
    "                 again. This double reading takes more time. If you want  ",
    "                 to avoid this, specify a maximum via this parameter.     ",
    "               =-1 Do not store records. Do not allocate memory for them. ",
    "                   The most complicated checking cannot be done. But all  ",
    "                   records are still read, converted to numbers, and      ",
    "                   checked individually as far as possible.               ",
    "               =-2 Do not store records. Do not allocate memory for them. ",
    "                   Do not convert to numbers. Do not check individually.  ",
    "                   This option can only be used with wtype=csvchop        ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "      scalco= check size of sx,sy,gx,gy after multiplying by this power   ",
    "              of 10 (1,10,100...). In this program all values from the    ",
    "              text file are checked to see if they fit in their SU keys.  ",
    "              The default in SUGEOMCSV is to multiply all coordinates     ",
    "              by 10 and so 10 is also the default here. But coordinates   ",
    "              are not actually multiplied by 10, instead the maximum size ",
    "              allowed is reduced by 10 during the checking performed here.",
    "              You only get warnings here, the numbers will still be       ",
    "              output to wfile.                                            ",
    "            * If you are confident your text files contain coordinates    ",
    "            * with only whole numbers, you can set this to 1 here and     ",
    "            * in SUGEOMCSV (but check first).                             ",
    "           ** This SEEMS like a problem about size, but it is actually a  ",
    "           ** problem about decimal digits. It you use 1 here you will    ",
    "           ** find that coordinates like 2223333.6 get rounded to 2223334 ",
    "           ** when SUGEOMCSV updates them to traces.                      ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "      scalel= check size of elevation and other related values after      ",
    "              multiplying by this power of 10 (1,10,100...). In this      ",
    "              program all values from the text file are checked to see if ",
    "              they fit in their SU keys. The default in SUGEOMCSV is to   ",
    "              multiply gelev,selev,sdepth,gdel,sdel,swdep,gwdep by 10 and ",
    "              so 10 is the default here. But these values are not actually",
    "              multiplied by 10, instead the maximum size allowed is       ",
    "              reduced by 10 during the checking performed herein.         ",
    "              You only get warnings here, the numbers will still be       ",
    "              output to wfile.                                            ",
    "            * If you are confident your text files contain these values   ",
    "            * with only whole numbers, you can set this to 1 here and     ",
    "            * in SUGEOMCSV (but check first).                             ",
    "           ** This SEEMS like a problem about size, but it is actually a  ",
    "           ** problem about decimal digits. It you use 1 here you will    ",
    "           ** find that values like 3333.6 get rounded to 3334            ",
    "           ** when SUGEOMCSV updates them to traces.                      ",
    "                                                                          ",
    "      outid=  Use this as first characters on output data records,        ",
    "              and on the output C_SU_SETID record.                        ",
    "              Default is the same value as the input setid.               ",
    "           =Q is the most likely id to change to because it is required   ",
    "              for (crooked) profiles input to suprofcsv and sunearcsv.    ",
    "                    Note: this value is automatically upper-cased unless  ",
    "                          you surround it by double-quotes.               ",
    "                          So q becomes Q unless you use double-quotes.    ",
    "                                                                          ",
    " ------------------------------------------------------------------------ ",
    " ------------------------------------------------------------------------ ",
    " ----------SPS rev2.1 Fixed Format Files---------------------------       ",
    "                                                                          ",
    "  For SPS format specification consult:                                   ",
    "    http://www.seg.org/resources/publications/misc/technical-standards    ",
    "                                                                          ",
    "  Remember that names= for fixed format files must be enclosed by leading ",
    "  and trailing numbers that specify their character ranges in the records.",
    "                                                                          ",
    "  The names=sps2 and forms=sps2 options expand as in the next examples.   ",
    "  For names=sps2all you get the same thing except that all NULL are       ",
    "  replaced with NUMB. And NUMB means the value is copied to output text   ",
    "  even though there is no SU name for it (or you will specify it later).  ",
    "  This means the value is put in the output file and will be available    ",
    "  when that file is input to a SpreadSheet or other program.              ",
    "                                                                          ",
    "     EXAMPLE for SPS rev2.1 X-file.                                       ",
    "     There are 15 fields defined in X-records but only 9 of them contain  ",
    "     values that have a reasonable chance of being useful (my opinion).   ",
    "     If you specify match=sps2 names=sps2 forms=sps2 setid=x              ",
    "     you get the following (as they appear in the output text file).      ",
    "       C_SU_MATCH,fldr,tracf                                              ",
    "       C_SU_SETID,X                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_ID,%.0f,%.0f,%.0f,%.0f,%.2f,%.2f,%.0f,,%.0f,%.0f,%.0f         ",
    "       C_SU_MORE,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.0f                       ",
    "       C_SU_NAMES                                                         ",
    "       C_SU_ID,2_null2_7,8_match1_15,16_null4_16,17_null5_17              ",
    "       C_SU_MORE,18_grnofr_27,28_grnlof_37,38_null8_38                    ",
    "       C_SU_MORE,39_matche1_cf_43,44_matche1_ct_48,49_matche1_ci_49       ",
    "       C_SU_MORE,50_grnors_59,60_gaps_rf_69,70_gaps_rt_79,80_null15_80    ",
    "                                                                          ",
    "     EXAMPLE for SPS rev2.1 S-file.                                       ",
    "     There are 18 fields defined in S-records but only 10 of them contain ",
    "     values that have a reasonable chance of being useful (my opinion).   ",
    "     If you specify match=sps2 names=sps2 forms=sps2 setid=s              ",
    "     you get the following (as they appear in the output text file).      ",
    "       C_SU_MATCH,grnofr,grnlof                                           ",
    "       C_SU_SETID,S                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_ID,%.2f,%.2f,%.0f,%.0f,%.0f,%.0f,%.1f,%.0f,%.0f,%.1f,%.1f     ",
    "       C_SU_MORE,%.1f,%.1f,%.0f,%.0f,%.0f,%.0f                            ",
    "       C_SU_NAMES                                                         ",
    "       C_SU_ID,2_grnofr_11,12_grnlof_21,22_null4_23,24_null5_24           ",
    "       C_SU_MORE,25_null6_26,27_sstat_30,31_sdepth_34,35_sdel_38,39_sut_40",
    "       C_SU_MORE,41_swdep_46,47_sx_55,56_sy_65,66_selev_71                ",
    "       C_SU_MORE,72_null15_74,75_null16_76,77_null17_78,79_null18_80      ",
    "                                                                          ",
    "     EXAMPLE for SPS rev2.1 R-file.                                       ",
    "     There are 18 fields defined in R-records but only 9 of them contain  ",
    "     values that have a reasonable chance of being useful (my opinion).   ",
    "     If you specify match=sps2 names=sps2 forms=sps2 setid=r              ",
    "     you get the following (as they appear in the output text file).      ",
    "       C_SU_MATCH,grnors,gaps                                             ",
    "       C_SU_SETID,R                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_ID,%.2f,%.2f,%.0f,%.0f,%.0f%.0f,%.1f,%.0f,%.0f,%.1f,%.1f      ",
    "       C_SU_MORE,%.1f,%.1f,%.0f,%.0f,%.0f,%.0f                            ",
    "       C_SU_NAMES                                                         ",
    "       C_SU_ID,2_grnors_11,12_gaps_21,22_null4_23,24_null5_24             ",
    "       C_SU_MORE,25_null6_26,27_gstat_30,31_null7_34,35_gdel_38,39_gut_40 ",
    "       C_SU_MORE,41_gwdep_46,47_gx_55,56_gy_65,66_gelev_71                ",
    "       C_SU_MORE,72_null15_74,75_null16_76,77_null17_78,79_null18_80      ",
    "                                                                          ",
    " Notes:                                                                   ",
    "        1. Since SPS records are only 80 characters wide, it is necessary ",
    "           to continue listing names on additional records (as above).    ",
    "           When needed, C_SU_MORE is put at start of continuation records.",
    "           It will be ignored (it does not count as a field).             ",
    "        2. Note C_SU_MORE is not permitted after C_SU_NAMES for           ",
    "           comma-separated files. (The actual data records will have a    ",
    "           certain number of commas, and the row containing the names     ",
    "           should/can/will have the same number of commas, so no need).   ",
    "        3. In the X-file specification above, I have used the SU keys     ",
    "           grnofr, grnlof, grnors, gaps to contain the values of          ",
    "           shot 3D line and point and receiver 3D line and point.         ",
    "           These are passed through the trace SU keys and used for        ",
    "           finding records within the S file and R file.                  ",
    "        4. Because grnofr, grnlof, grnors, gaps are short integers        ",
    "           you can only have line and point number magnitudes less than   ",
    "           roughtly 32765.                                                ",
    "       *** See the process= options for what you can do if you encounter  ",
    "           this numeric magnitude limitation.                             ",
    "                                                                          ",
    "                                                                          ",
    " c_su_id IS SPECIAL ********************                                  ",
    "                                                                          ",
    "  (a) Its character range is taken from the length of the value specified ",
    "      for setid (for instance S has range 1 to 1, FRED has range 1 to 4). ",
    "  (b) The value specified for id is case-sensitive (r is not R). The id   ",
    "      value is the only thing case-sensitive in this program except for   ",
    "      the file names. So this program does not care if you use parameter  ",
    "      records starting with C_SU or c_su, but you want other programs to  ",
    "      ignore C_SU records, so a capital C is better.                      ",
    "                                                                          ",
    "                                                                          ",
    " *** Special names ***                                                    ",
    "                                                                          ",
    " You may have noticed the names match1 and matche1 in X-file names= list. ",
    " This means substitute this with the corresponding name from match= list. ",
    " Where match1 means the first match= name and match2 means second match=  ",
    " and so on. And where matche1 means ending match= name and matche2 means  ",
    " the second-to-ending match= name and so on.                              ",
    " This facility exists because values like field record numbers do not     ",
    " always end up in the fldr SU key in the input traces. Similarly, values  ",
    " like channel number do not always end up in the tracf SU key. These      ",
    " substitutions allow you to specify a more generic setup.                 ",
    " (They can only have a single digit from 1 to 9 on the end).              ",
    "                                                                          ",
    " Further, note the extra _cf _ct _ci and _rf _rt on the ends of some of   ",
    " the X-file names list. The first 3 indicate that these are the from, to, ",
    " and increment for the channel ranges and the second 2 are the from, to   ",
    " for the corresponding receiver range. These extras are how this program  ",
    " recognizes which values to use to compute the output receiver value from ",
    " the input channel value (again, this is generic, this program does not   ",
    " care if these are actual channels and receivers).                        ",
    "                                                                          ",
    " ------------------------------------------------------------------------ ",
    " ----------SPS rev1 Fixed Format Files-----------------------------       ",
    "                                                                          ",
    " I include the SPS rev1 format example because you are quite likely to    ",
    " encounter it. It also illustrates an important issue. The rev1 format    ",
    " allows the line values to be alphanumeric. Therefore 14_grnofr_29 and    ",
    " and 48_grnors_63 specified below are likely to cause error-warnings in   ",
    " this program because it will attempt to read these fields as numbers.    ",
    "  - You can use the process= options (such as trimz) to alter the file.   ",
    "  - Or, you can change the range to only include the numeric part of      ",
    "    the line name (example: change 14_grnofr_29 to 20_grnofr_24).         ",
    "  - Or, for 2D lines, you do not actually need the line names so you      ",
    "    can change 14_grnofr_29 to 14_null_29.                                ",
    " But all changes will need to be done for the corresponding values in the ",
    " S and R files (you usually encounter these problems in the X file first).",
    "                                                                          ",
    "     EXAMPLE for SPS rev1 X-file.                                         ",
    "       C_SU_MATCH,fldr,tracf                                              ",
    "       C_SU_SETID,X                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_%s,%.0f,%.0f,%.0f,%.0f,%.2f,%.2f,%.0f,,%.0f,%.0f,%.0f         ",
    "       C_SU_MORE,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.0f                       ",
    "       C_SU_NAMES                                                         ",
    "       C_SU_ID,2_null_7,8_match1_11,12_null_12,13_null_13                 ",
    "       C_SU_MORE,14_grnofr_29,30_grnlof_37,38_null_38                     ",
    "       C_SU_MORE,39_matche1_cf_42,43_matche1_ct_46,47_matche1_ci_47       ",
    "       C_SU_MORE,48_grnors_63,64_gaps_rf_71,72_gaps_rt_79,80_null_80      ",
    "                                                                          ",
    "     EXAMPLE for SPS rev1 S-file.                                         ",
    "       C_SU_MATCH,grnofr,grnlof                                           ",
    "       C_SU_SETID,S                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_%s,%.2f,%.2f,%.0f,%.0f,%.0f,%.0f,%.1f,%.0f,%.0f,%.1f,%.1f     ",
    "       C_SU_MORE,%.1f,%.1f,%.0f,%.0f,%.0f,%.0f                            ",
    "       C_SU_NAMES                                                         ",
    "       C_SU_ID,2_grnofr_17,18_grnlof_25,26_null_26,26_null_26             ",
    "       C_SU_MORE,27_null_28,29_sstat_32,33_sdepth_36,37_sdel_40,41_sut_42 ",
    "       C_SU_MORE,43_swdep_46,47_sx_55,56_sy_65,66_selev_71                ",
    "       C_SU_MORE,72_null_74,75_null_76,77_null_78,79_null_80              ",
    "                                                                          ",
    "     EXAMPLE for SPS rev1 R-file.                                         ",
    "       C_SU_MATCH,grnors,gaps                                             ",
    "       C_SU_SETID,R                                                       ",
    "       C_SU_FORMS                                                         ",
    "       C_SU_%s,%.2f,%.2f,%.0f,%.0f,%.0f%.0f,%.1f,%.0f,%.0f,%.1f,%.1f      ",
    "       C_SU_MORE,%.1f,%.1f,%.0f,%.0f,%.0f,%.0f                            ",
    "       C_SU_NAMES                                                         ",
    "       C_SU_ID,2_grnors_17,18_gaps_25,26_null_26,26_null_26               ",
    "       C_SU_MORE,27_null_28,29_gstat_32,33_null_36,37_gdel_40,41_gut_42   ",
    "       C_SU_MORE,43_gwdep_46,47_gx_55,56_gy_65,66_gelev_71                ",
    "       C_SU_MORE,72_null_74,75_null_76,77_null_78,79_null_80              ",
    "                                                                          ",
    "   When parsing fixed format files it is not actually needed to specify   ",
    "   fields that you decide to null. But it is easier to make changes in    ",
    "   the future if you retain all defined fields in the names= lists.       ",
    "   That is, if you do not want 29_gstat_32 from the R-file, I advise      ",
    "   changing it to 29_null_32 rather than removing it from the list.       ",
    "                                                                          ",
    "   The keen observer will notice that 26_null_26 appears twice in the     ",
    "   SPS1 S and R names lists. This more-easily allows conversion from      ",
    "   SPS1 to SPS2 (and vice-versa).  Example usage to convert SPS1 to SPS2: ",
    "       sutoolcsv rfile=A.txt wfile=B.csv setid=x                          ",
    "                 match=sps1 names=sps1 forms=sps1                         ",
    "       sutoolcsv rfile=B.csv wfile=C.txt setid=x                          ",
    "                 names=sps2 forms=sps2 forms=sps2                         ",
    "   This works because there are the same number of fields (commas) made   ",
    "   by options SPS1 and SPS2 and the values are in the same order.         ",
    "   *** But going from sps1 to sps2 it is almost certain you will have to  ",
    "   *** modify A.txt before getting good results. See process= options.    ",
    "                                                                          ",
    "                                                                          ",
    "                                                                          ",
    " ---Special Consideration For Relational Files (such as SPS X-files)------",
    "                                                                          ",
    "    Part of using Relational files (like SPS X-file) is similar to using  ",
    "    simpler text files. Relational records have values which match one or ",
    "    more key values from the input trace header. Usually this is something",
    "    called the field record number (fldr). Its values from the header are ",
    "    then searched for in the X-records using the field with the same name.",
    "    But usually multiple X-records match each fldr number. For instance,  ",
    "    for a split-spread 2D, you should expect 2 X-records for each fldr.   ",
    "    These 2 X-records describe which channels are on which receiver points",
    "    for each particular shot. For 3D surveys, there are often 10 or so    ",
    "    X-records for each fldr. For example, a fldr may record 2400 traces   ",
    "    with the recording geophone layout pattern extending over 10 receiver ",
    "    lines with 240 receiver points each. The next fldr will usually also  ",
    "    record 10 lines by 240 points each, but they may be a different set   ",
    "    of 10 lines and a differnt set of 240 points. So, SUGEOMCSV program   ",
    "    must read a fldr AND a channel number from the input trace header.    ",
    "    The channel number allows the SUGEOMCSV program to determine which of ",
    "    the 10 X-records to use for that input trace because each of the      ",
    "    10 X-records contains a channel range (from,to,increment). Each of    ",
    "    the 10 X-records also has a receiver point range (from,to). These two ",
    "    ranges allow computation of the specific receiver point number        ",
    "    corresponding to the channel number from the trace header.            ",
    "    The from,to,increment channel names are seen in the X-record examples ",
    "    above (39_matche1_cf_43,44_matche1_ct_48,49_matche1_ci_49) and need to",
    "    have _cf _ct _ci to identify them. Similarly, the from,to receiver    ",
    "    point names are seen in example above (60_gaps_rf_69,70_gaps_rt_79).  ",
    "    Notice also the names match1 and matche1 in the X-file name list. This",
    "    means substitute this with the corrasponding name from the match= list",
    "    (with matche1 meaning the name on the end of the match list).         ",
    "    This facility exists because values like field record numbers do not  ",
    "    always end up in the fldr SU key in input trace headers. Similarly,   ",
    "    values like channel number do not always end up in the tracf SU key.  ",
    "    Substitutions allow you to specify a more generic setup for X-files   ",
    "    (but this is general, match= substitution is not limited to X-files). ",
    "                                                                          ",
    "    If the names= list does not indicate this is a Relational file then   ",
    "    it is assumed to be similar to S and R files. So this program         ",
    "    presumes that a match= list with 2 names are lines and points on those",
    "    lines. It then warns about inconsistant incrementing of those values. ",
    "    That is, it is trying to warn you if some point records are missing   ",
    "    within the lines. If there is only 1 name on the match= list then this",
    "    program presumes those are points. If the input file is not similar   ",
    "    to S and R files, you will get some false warnings.                   ",
    "                                                                          ",
    " -----------------------------------------------------------------        ",
    "                                                                          ",
];

#[derive(Clone)]
pub(crate) struct PointInfo {
    dfield: Vec<f64>,
    lfield: Vec<i64>,
}

// Borrow the owned fields using the same two-pointer element size as the C
// PointInfo. qsort may permute equal keys, and the layout checks observe that
// order. The owning PointInfo array remains alive and unchanged during sorting
// and checking; only these Copy views move.
#[repr(C)]
#[derive(Clone, Copy)]
pub(crate) struct PointView<'a> {
    dfield: &'a Vec<f64>,
    lfield: &'a Vec<i64>,
}

unsafe extern "C" {
    fn qsort(
        base: *mut std::ffi::c_void,
        count: usize,
        size: usize,
        compare: unsafe extern "C" fn(*const std::ffi::c_void, *const std::ffi::c_void) -> i32,
    );
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    if !isatty(STDIN) || !isatty(STDOUT) {
        err!("**** Error: Traces are not input and not output by this program.");
    }

    let mut nicerecord = 1i32;
    getparint("nicerecord", &mut nicerecord);
    if nicerecord < 1 {
        err!("**** Error: nicerecord= cannot be less than 1");
    }
    let mut maxrecords = 0i32;
    getparint("maxrecords", &mut maxrecords);
    if maxrecords < -2 {
        err!("**** Error: maxrecords= cannot be less than -2");
    }
    let mut ftol = 0.01f32;
    getparfloat("tolr", &mut ftol);
    if (ftol as f64) < 0.000000001 {
        err!("**** Error: tolr= must be larger.");
    }
    let mut dtol = ftol as f64;
    let dtolh = dtol / 2.;
    dtol = 1.0 / dtol;
    let mut unrepeat = -2147483645i32;
    getparint("unrepeat", &mut unrepeat);
    if maxrecords == -2 && unrepeat > -2147483645 {
        err!("**** Error: maxrecords=-2 and unrepeat not allowed at same time.");
    }
    let mut Rname = String::new();
    if getparstring("rfile", &mut Rname) == 0 {
        err!("**** Error: rfile= text file name must be specified.");
    }
    let mut rdel = b',';
    let mut tdel = String::new();
    let has_tdel = getparstring("rdelim", &mut tdel) != 0;
    if has_tdel {
        if tdel.len() == 1 {
            rdel = tdel.as_bytes()[0];
        } else if tdel.len() == 3 && tdel.starts_with('"') && tdel.ends_with('"') {
            rdel = tdel.as_bytes()[1];
        } else {
            err!("**** Error: rdelim= specification not recognized.");
        }
    }
    let mut irtype = 1;
    let mut rtype = String::new();
    if countparval("rtype") > 0 {
        getparstring("rtype", &mut rtype);
        rtype.make_ascii_lowercase();
        if rtype == "csv" {
            irtype = 1;
        } else if rtype == "fixed" {
            irtype = 0;
        } else {
            err!("**** Error: rtype= option not recognized.");
        }
    } else if Rname.len() > 2 {
        irtype = if Rname.ends_with("csv") { 1 } else { 0 };
    }
    if irtype == 1 {
        if !has_tdel {
            rdel = b',';
        }
    } else if has_tdel {
        err!("**** Error: you cannot specify rdelim= if rtype=fixed.");
    }
    let fpR = fopen(&Rname, "r").unwrap_or_else(|| {
        err!("**** Error opening the rfile text file.");
    });
    let mut match_ = Vec::<String>::new();
    let match_count = countparval("match") as usize;
    if match_count > 0 {
        getparstringarray("match", &mut match_);
    }
    let mut Rid = String::new();
    let mut lenid = 0usize;
    if countparval("setid") > 0 {
        getparstring("setid", &mut Rid);
        lenid = Rid.len();
    }
    let mut names = Vec::<String>::new();
    let names_count = countparval("names") as usize;
    if names_count > 0 {
        getparstringarray("names", &mut names);
        if names.len() != names_count {
            vendor_segv();
        }
        for name in &mut names {
            name.make_ascii_lowercase();
        }
    }
    let mut forms = Vec::<String>::new();
    let forms_count = countparval("forms") as usize;
    if forms_count > 0 {
        getparstringarray("forms", &mut forms);
    }
    let mut Wname = String::new();
    if getparstring("wfile", &mut Wname) == 0 {
        err!("**** Error: wfile= output text file name must be specified.");
    }
    if Rname == Wname {
        err!("**** Error: wfile= output file must be different than input.");
    }
    let mut iwtype = 1;
    let mut iwchop = 0;
    let mut wtype = String::new();
    if countparval("wtype") > 0 {
        getparstring("wtype", &mut wtype);
        wtype.make_ascii_lowercase();
        if wtype == "csv" {
            iwtype = 1;
        } else if wtype == "csvchop" {
            iwtype = 1;
            iwchop = 1;
        } else if wtype == "fixed" {
            iwtype = 0;
        } else {
            err!("**** Error: wtype= option not recognized.");
        }
    } else if Wname.len() > 2 {
        iwtype = if Wname.ends_with("csv") { 1 } else { 0 };
    }
    if maxrecords == -2 && iwchop == 0 {
        err!("**** Error: maxrecords=-2 only allowed for wtype=csvchop.");
    }
    let mut iwidth = 0i32;
    let mut iC_SU = 1;
    getparint("width", &mut iwidth);
    if iwidth > 0 {
        if iwtype != 0 {
            err!("**** Error: positive width= only allowed for wtype=fixed.");
        }
        if iwidth < 31 {
            err!("**** Error: positive width= has to be greater or equal to 30.");
        }
    }
    if iwidth < 0 {
        if iwidth == -1 {
            iwidth = 0;
        } else {
            if iwtype != 0 {
                err!("**** Error: only width=-1 allowed when wtype is not fixed.");
            }
            iwidth = iwidth.wrapping_neg();
        }
        iC_SU = 0;
    }
    let mut pross = Vec::<String>::new();
    let pross_count = countparval("process") as usize;
    if pross_count > 0 {
        if irtype != 0 {
            err!("**** Error: you can only specify process= when rfile type is fixed.");
        }
        getparstringarray("process", &mut pross);
        if pross.len() != pross_count {
            vendor_segv();
        }
        for s in &mut pross {
            s.make_ascii_lowercase();
        }
    }
    let mut iscalel = 10i32;
    getparint("scalel", &mut iscalel);
    if iscalel == -1 {
        iscalel = 1;
    }
    let nscalel = iscalel.wrapping_abs();
    if ![1, 10, 100, 1000, 10000, 100000, 1000000, 10000000].contains(&nscalel) {
        err!("**** Error: scalel= must be signed powers of 10 (1,10,100...-10,-100,...)");
    }
    let dscalel = if iscalel > 0 {
        iscalel as f64
    } else {
        -1. / iscalel as f64
    };
    let mut iscalco = 10i32;
    getparint("scalco", &mut iscalco);
    if iscalco == -1 {
        iscalco = 1;
    }
    let nscalco = iscalco.wrapping_abs();
    if ![1, 10, 100, 1000, 10000, 100000, 1000000, 10000000].contains(&nscalco) {
        err!("**** Error: scalco= must be signed powers of 10 (1,10,100...-10,-100,...)");
    }
    let dscalco = if iscalco > 0 {
        iscalco as f64
    } else {
        -1. / iscalco as f64
    };

    let maxtext = 10001usize;
    let mut textraw = vec![0u8; maxtext];
    let mut textbeg = vec![0u8; maxtext];
    let mut textraw2 = vec![0u8; maxtext];
    let mut names_more = 0;
    let mut forms_more = 0;
    if names_count == 0 || forms_count == 0 || match_count == 0 || lenid == 0 {
        let in_num_names = names_count;
        let in_num_forms = forms_count;
        let in_num_to_sort_by = match_count;
        let in_lenid = lenid;
        let mut num_c_su_names = 0;
        let mut num_c_su_forms = 0;
        let mut num_c_su_match = 0;
        let mut num_c_su_setid = 0;
        let mut read_names = 0;
        let mut read_forms = 0;
        loop {
            // Read the C record limit in bytes. Core fgets returns a lossy
            // UTF-8 String, which can expand beyond this fixed-size buffer.
            let mut nread = 0;
            while nread < maxtext - 1 {
                let c = fgetc(&fpR);
                if c == EOF {
                    break;
                }
                textraw[nread] = c as u8;
                nread += 1;
                if c == b'\n' as i32 {
                    break;
                }
            }
            if nread == 0 {
                break;
            }
            textraw[nread] = 0;
            if (in_num_names > 0 || read_names == -2)
                && (in_num_forms > 0 || read_forms == -2)
                && (match_count > 0 || !match_.is_empty())
                && lenid > 0
            {
                break;
            }
            let mut tsize = 0;
            for &b in &textraw[..nread] {
                if matches!(b, 0 | b'\n' | b'\r') {
                    break;
                }
                if b != b' ' && b != b'\t' {
                    textbeg[tsize] = b;
                    tsize += 1;
                }
            }
            if lenid == 0 {
                textbeg[..10].make_ascii_lowercase();
                if textbeg.starts_with(b"c_su_setid") {
                    textbeg[tsize] = 0;
                    let ids = tparse(&textbeg[..tsize], b',');
                    Rid = ids.get(1).cloned().unwrap_or_default();
                    lenid = Rid.len();
                }
            }
            // sizeof(textbeg) in C is the full array; keep stale suffix bytes.
            textbeg.make_ascii_lowercase();
            if textbeg.starts_with(b"c_su_match") {
                num_c_su_match += 1;
            }
            if textbeg.starts_with(b"c_su_setid") {
                num_c_su_setid += 1;
            }
            if textbeg.starts_with(b"c_su_names") {
                num_c_su_names += 1;
            }
            if textbeg.starts_with(b"c_su_forms") {
                num_c_su_forms += 1;
            }
            if match_count == 0 && match_.is_empty() && textbeg.starts_with(b"c_su_matc") {
                textbeg[tsize] = 0;
                let found = tparse(&textbeg[..tsize], b',');
                for s in found.into_iter().skip(1) {
                    if s.starts_with("null") {
                        break;
                    }
                    match_.push(s);
                }
            }
            if read_names == -1 {
                if textbeg.starts_with(b"c_su_more") {
                    read_names = 1;
                    names_more = 1;
                } else {
                    read_names = -2;
                }
            }
            if read_names > 0 {
                textbeg[tsize] = 0;
                names.extend(tparse(&textbeg[..tsize], b','));
                read_names = -1;
            }
            if in_num_names == 0 && textbeg.starts_with(b"c_su_names") {
                read_names = 1;
            }
            if read_forms == -1 {
                if textbeg.starts_with(b"c_su_more") {
                    read_forms = 1;
                    forms_more = 1;
                } else {
                    read_forms = -2;
                }
            }
            if read_forms > 0 {
                textbeg[tsize] = 0;
                forms.extend(tparse(&textbeg[..tsize], b','));
                read_forms = -1;
            }
            if in_num_forms == 0 && textbeg.starts_with(b"c_su_forms") {
                read_forms = 1;
            }
        }
        fseeko(&fpR, 0, SEEK_SET);
        let mut ierr = 0;
        if in_num_names < 1 {
            if num_c_su_names > 1 {
                ierr = 1;
                warn!("**** Error: Text file has more than one C_SU_NAMES parameter record.");
                warn!("****        Remove duplicates or override with names= on command line.");
            } else if num_c_su_names == 0 {
                ierr = 1;
                warn!("**** Error: Text file has no C_SU_NAMES parameter record.");
                warn!("****        Add it, or specify names= on command line.");
            }
        }
        if in_num_forms < 1 {
            if num_c_su_forms > 1 {
                ierr = 1;
                warn!("**** Error: Text file has more than one C_SU_FORMS parameter record.");
                warn!("****        Remove duplicates or override with forms= on command line.");
            } else if num_c_su_forms == 0 {
                ierr = 1;
                warn!("**** Error: Text file has no C_SU_FORMS parameter record.");
                warn!("****        Add it, or specify forms= on command line.");
            }
        }
        if in_num_to_sort_by < 1 {
            if num_c_su_match > 1 {
                ierr = 1;
                warn!("**** Error: Text file has more than one C_SU_MATCH record.");
                warn!("****        Remove duplicates or override with match= on command line.");
            } else if num_c_su_match == 0 {
                ierr = 1;
                warn!("**** Error: Text file has no C_SU_MATCH record.");
                warn!("****        Add it, or specify match= command line.");
            }
        }
        if in_lenid == 0 {
            if num_c_su_setid > 1 {
                ierr = 1;
                warn!("**** Error: Text file has more than one C_SU_SETID record.");
                warn!("****        Remove duplicates or override with setid= on command line.");
            } else if num_c_su_setid == 0 {
                ierr = 1;
                warn!("**** Error: Text file has no C_SU_SETID record.");
                warn!("****        Add it, or specify setid= on command line.");
            }
        }
        if ierr > 0 {
            err!("**** Error: Text file has duplicate or missing C_SU_NAMES or C_SU_MATCH or C_SU_SETID records.");
        }
    }
    // Resolve setid and outid. Quotes embedded in the argument preserve case.
    if lenid > 3 && Rid.starts_with('"') && Rid.ends_with('"') {
        // C shifts the opening quote away without moving the terminator.
        // Keep the two closing quotes in storage; lenid limits comparisons.
        Rid.remove(0);
        Rid.push('"');
        lenid -= 2;
    } else {
        Rid.make_ascii_uppercase();
    }
    let mut isetid = 1;
    if lenid == 3 && Rid[..lenid].eq_ignore_ascii_case("any") {
        lenid = 0;
    } else if lenid == 4 && Rid[..lenid].eq_ignore_ascii_case("none") {
        lenid = 0;
        isetid = 0;
    }
    let mut lenod = lenid;
    let mut Rod = Rid.clone();
    if countparval("outid") > 0 {
        getparstring("outid", &mut Rod);
        lenod = Rod.len();
        if lenod > 3 && Rod.starts_with('"') && Rod.ends_with('"') {
            Rod.remove(0);
            Rod.push('"');
            lenod -= 2;
        } else {
            Rod.make_ascii_uppercase();
        }
    }

    if match_count > 0 && match_.len() != match_count {
        vendor_segv();
    }
    if match_.first().is_some_and(|s| s == "sps2" || s == "sps1") {
        match Rid.as_str() {
            "X" => match_ = vec!["fldr".into(), "tracf".into()],
            "S" => match_ = vec!["grnofr".into(), "grnlof".into()],
            "R" => match_ = vec!["grnors".into(), "gaps".into()],
            _ => {}
        }
    }
    let num_to_sort_by = match_.len();
    let mut namesps = 0;
    if names.first().is_some_and(|s| s == "sps2" || s == "sps2all") {
        namesps = 2;
        if let Some(a) = presets::names(2, &Rid, names[0] == "sps2all") {
            names = a;
        }
    }
    if names.first().is_some_and(|s| s == "sps1" || s == "sps1all") {
        namesps = 1;
        if let Some(a) = presets::names(1, &Rid, names[0] == "sps1all") {
            names = a;
        }
    }
    if namesps != 0 && Rid == "X" {
        // In the vendor build the three 16-character names are allocated
        // in adjacent 16-byte slots. Each following strcpy overwrites the
        // preceding terminator, joining fields 8 through 11. Preserve the
        // resulting deterministic parse error, including its full text.
        names[8] = names[8..=11].concat();
    }
    let mut ilead = vec![0i32; names.len()];
    let mut itrail = vec![0i32; names.len()];
    let mut ncase = vec![0i32; names.len()];
    let mut nspot = vec![0usize; names.len()];
    let mut valmx = vec![0.0f64; names.len()];
    let mut kcase = vec![0i32; num_to_sort_by];
    let mut ktol = vec![0usize; num_to_sort_by];
    let mut formsps = 0;
    if forms_count > 0 && forms.len() != forms_count {
        vendor_segv();
    }
    if forms.first().is_some_and(|s| s == "sps2" || s == "sps1") {
        formsps = if forms[0] == "sps1" { 1 } else { 2 };
        if let Some(a) = presets::forms(&Rid) {
            forms = a;
        }
    }
    if namesps != formsps {
        warn!("Warning: Different sps options for names= and forms=. Unusual, but sometimes intentional.");
    }
    forms.retain(|s| !s.starts_with("c_su_more"));
    let mut form2 = forms.clone();
    names.retain(|s| !s.starts_with("c_su_more"));
    let mut name2 = names.clone();
    let num_names = names.len();
    let num_forms = forms.len();
    for _ in num_forms..num_names {
        forms.push(forms[num_forms - 1].clone());
        form2.push(forms[num_forms - 1].clone());
    }
    if names.first().is_some_and(|s| s.starts_with("c_su_id")) {
        names[0] = "null".into();
    }
    let mut namex = vec![String::new(); num_names];
    let mut incomma = 0;
    let mut extra_parts = 0;
    let mut maxtrail = 0;
    for n in 0..num_names {
        let nparts = tparse(names[n].as_bytes(), b'_');
        ilead[n] = -1;
        itrail[n] = -1;
        match nparts.len() {
            1 => names[n] = nparts[0].clone(),
            2 => {
                names[n] = nparts[0].clone();
                namex[n] = nparts[1].clone();
            }
            3 | 4 => {
                ilead[n] = atoi(&nparts[0]);
                itrail[n] = atoi(&nparts[nparts.len() - 1]);
                if itrail[n] < ilead[n] || ilead[n] < 1 {
                    err!(
                        "**** Error: Your names= list has an entry with incorrect integers: %s",
                        names[n]
                    );
                }
                if itrail[n] > maxtrail {
                    maxtrail = itrail[n];
                }
                names[n] = nparts[1].clone();
                if nparts.len() == 4 {
                    namex[n] = nparts[2].clone();
                }
            }
            _ => err!(
                "**** Error: Your names= list has an entry that parses incorrectly: %s",
                names[n]
            ),
        }
        if !names[n].starts_with("null") && names[n] != "c_su_id" {
            if nparts.len() < 3 {
                incomma = 1;
            }
            if nparts.len() == 2 || nparts.len() == 4 {
                extra_parts += 1;
            }
        }
    }
    let mut lerr = 0;
    if incomma == 1 && irtype == 0 {
        lerr = 1;
        warn!("The rtype=fixed but at least one non-null name has no leading and trailing range.");
    }
    if extra_parts != 0 && extra_parts != 5 {
        lerr = 1;
        warn!("**** Error: Your names= list only has some of _cf _ct _ci _rf _rt. Need all 5 or none.");
    }
    if incomma == 1 && names_more == 1 {
        lerr = 1;
        warn!("**** Error: C_SU_MORE not permitted after C_SU_NAMES for comma-separated files.");
    }
    if incomma == 1 && forms_more == 1 {
        lerr = 1;
        warn!("**** Error: C_SU_MORE not permitted after C_SU_FORMS for comma-separated files.");
    }
    for n in 0..num_names {
        if !names[n].starts_with("null") {
            for k in 1..=9 {
                if names[n] == format!("match{k}") && num_to_sort_by >= k {
                    names[n] = match_[k - 1].clone();
                }
                if names[n] == format!("matche{k}") && num_to_sort_by >= k {
                    names[n] = match_[num_to_sort_by - k].clone();
                }
            }
            if names[n].starts_with("match") {
                lerr = 1;
                warn!(
                    "**** Error: Name  %s  could not be substituted from match= list.",
                    names[n]
                );
            }
            for m in n + 1..num_names {
                if names[n] == names[m] && namex[n] == namex[m] {
                    lerr = 1;
                    warn!(
                        "**** Error: Name  %s  exists at least twice in the names list.",
                        names[n]
                    );
                }
            }
        }
    }
    if lerr > 0 {
        err!("**** Error: Related to names= or C_SU_NAMES record (details above).");
    }
    if iwidth == 0 {
        iwidth = maxtrail;
        if iwidth < 80 {
            iwidth = 80;
        }
    }
    if iwidth < maxtrail && iwtype == 0 {
        err!("**** Error: Your specified width= not wide enough for maximum trailing integer on a name.");
    }
    for i in 0..num_to_sort_by {
        kcase[i] = GetCase(&match_[i]);
        if kcase[i] < 1 {
            err!(
                "**** Error: a match name not recognized (or not allowed). %s",
                match_[i]
            );
        }
    }

    let mut plead = vec![-1i32; pross.len()];
    let mut ptrail = vec![-1i32; pross.len()];
    let mut pflag = vec![0i32; pross.len()];
    let mut pvalu = vec![0.0f64; pross.len()];
    for n in 0..pross.len() {
        let nparts = tparse(pross[n].as_bytes(), b'_');
        pross[n].clear();
        if nparts.len() == 3 || nparts.len() == 4 {
            plead[n] = atoi(&nparts[0]);
            ptrail[n] = atoi(&nparts[nparts.len() - 1]);
            if ptrail[n] < plead[n] || plead[n] < 1 {
                err!(
                    "**** Error: Your process= list has an entry with incorrect integers: %s",
                    pross[n]
                );
            }
            pross[n] = nparts[1].clone();
            if nparts.len() == 4 {
                let (v, used) = strtod(&nparts[2]);
                if used > 0 {
                    pvalu[n] = v;
                }
            }
        } else {
            err!(
                "**** Error: Your process= list has some entries that do not parse understandably."
            );
        }
        if pross[n].starts_with("trimz") {
            pflag[n] = 3;
        } else if pross[n].starts_with("trim") {
            pflag[n] = 2;
        } else if pross[n].starts_with("sub") {
            pflag[n] = -1;
            pvalu[n] *= -1.;
        } else if pross[n].starts_with("add") {
            pflag[n] = -1;
        } else if pross[n].starts_with("zero") {
            pflag[n] = 0;
        } else if pross[n].starts_with("blank") {
            pflag[n] = 1;
        } else if pross[n].starts_with("div") {
            pflag[n] = -2;
            if pvalu[n] == 0.0 {
                err!("**** Error: Your process= list says to divide by 0.");
            }
            pvalu[n] = 1.0 / pvalu[n];
        } else if pross[n].starts_with("mul") {
            pflag[n] = -2;
        } else {
            err!("**** Error: Your process= list has an option that is not recognized.");
        }
    }
    let mut numcases = 0usize;
    for i in 0..num_names {
        let iamcase = GetCase(&names[i]);
        if iamcase < 0 {
            err!(
                "**** Error: Name  %s  in the names list is not recognized.",
                names[i]
            );
        } else if iamcase > 0 {
            ncase[numcases] = iamcase;
            ilead[numcases] = ilead[i];
            itrail[numcases] = itrail[i];
            names[numcases] = names[i].clone();
            if iwtype != 0 {
                name2[numcases] = name2[i].clone();
                form2[numcases] = form2[i].clone();
            }
            namex[numcases] = namex[i].clone();
            forms[numcases] = forms[i].clone();
            nspot[numcases] = if irtype == 1 { i } else { numcases };
            if names[i].starts_with("numb") {
                valmx[numcases] = 1.7e307;
            } else {
                valmx[numcases] = match hdtype(&names[i]) {
                    "i" => 2147483645.,
                    "h" => 32765.,
                    "u" => 65533.,
                    "f" => 3.4e37,
                    _ => 1.7e307,
                };
            }
            if ["gelev", "selev", "sdepth", "gdel", "sdel", "swdep", "gwdep"]
                .contains(&names[i].as_str())
            {
                valmx[numcases] /= dscalel;
            }
            if ["sx", "sy", "gx", "gy"].contains(&names[i].as_str()) {
                valmx[numcases] /= dscalco;
            }
            numcases += 1;
        }
    }
    let mut mapx = [-1i32; 10];
    let mut kerr = 0;
    if extra_parts == 5 {
        for i in 0..numcases {
            let slot = match namex[i].as_str() {
                "cf" => Some(0),
                "ct" => Some(1),
                "ci" => Some(2),
                "rf" => Some(3),
                "rt" => Some(4),
                _ => None,
            };
            if let Some(slot) = slot {
                if mapx[slot] != -1 {
                    kerr = 1;
                    match slot {
                        0 => warn!("**** Error: Only one name with _cf appended is allowed."),
                        1 => warn!("**** Error: Only one name with _ct appended is allowed."),
                        2 => warn!("**** Error: Only one name with _ci appended is allowed."),
                        3 => warn!("**** Error: Only one name with _rf appended is allowed."),
                        _ => warn!("**** Error: Only one name with _rt appended is allowed."),
                    }
                }
                mapx[slot] = i as i32;
                if slot == 1 || slot == 2 || slot == 4 {
                    ncase[i] = 0;
                }
            } else if names[i] == "grnofr" {
                mapx[7] = i as i32;
            } else if names[i] == "grnlof" {
                mapx[8] = i as i32;
            }
        }
        if mapx[..5].iter().all(|i| *i >= 0) {
            if names[mapx[0] as usize] != names[mapx[1] as usize]
                || names[mapx[0] as usize] != names[mapx[2] as usize]
            {
                kerr = 1;
                warn!("**** Error: _cf _ct _ci must be appended to the same name.");
            }
            if names[mapx[3] as usize] != names[mapx[4] as usize] {
                kerr = 1;
                warn!("**** Error: _rf _rt must be appended to the same name.");
            }
        } else {
            // C indexes before the names array when appendices are missing.
            kerr = 1;
        }
        if kerr > 0 {
            err!("**** Error: Your _cf _ct _ci _rf _rt specification is incorrect.");
        }
        for k in 0..num_to_sort_by {
            if names[mapx[0] as usize] == match_[k] {
                mapx[5] = k as i32;
                mapx[6] = kcase[k];
            }
        }
        if mapx[5] == -1 {
            err!("**** Error: Cannot find the name appended with _cf _ct and _ci in match= list.");
        }
    }
    for k in 0..num_to_sort_by {
        let mut ifound = 0;
        for n in 0..numcases {
            if ncase[n] == kcase[k] {
                ktol[k] = n;
                ifound = 1;
                break;
            }
        }
        if ifound == 0 {
            err!("**** Error: match= name not found in names=");
        }
    }
    for n in 0..numcases {
        for k in 0..num_to_sort_by {
            if ncase[n] == kcase[k] {
                ncase[n] = 0;
            }
        }
    }
    let numR: usize;
    if maxrecords == 0 {
        numR = countRec(&fpR, &mut textraw, &Rid, lenid, nicerecord);
        warn!("Counted %d data records.", numR as i32);
        if numR == 0 {
            err!("**** No data records found. Wrong setid value? Wrong file?");
        }
    } else if maxrecords > 0 {
        numR = maxrecords as usize;
    } else {
        numR = 1;
    }
    let mut RecInfo = vec![
        PointInfo {
            dfield: vec![0.0; numcases],
            lfield: vec![0; num_to_sort_by]
        };
        numR
    ];
    let fpW = fopen(&Wname, "w").unwrap_or_else(|| {
        err!("**** Error opening the wfile output text file.");
    });
    textraw.fill(0);
    textraw2.fill(0);
    textbeg.fill(0);

    let iwidth = iwidth as usize;
    if iC_SU > 0 {
        if iwtype == 0 {
            textraw[..iwidth].fill(b' ');
            textraw[iwidth] = b'\n';
        }
        textraw[..11].copy_from_slice(b"C_SU_MATCH\0");
        let mut mspot = 9usize;
        for key in &match_ {
            mspot += 1;
            textraw[mspot] = b',';
            let mleng = key.len();
            textraw[mspot + 1..mspot + 1 + mleng].copy_from_slice(key.as_bytes());
            mspot += mleng;
        }
        if iwtype == 0 && mspot > iwidth {
            err!("**** Error: too many match= for width= size.");
        }
        if iwtype != 0 {
            textraw[mspot + 1] = b'\n';
            textraw[mspot + 2] = 0;
        }
        let end = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
        fputs_bytes(&textraw[..end], &fpW);

        if iwtype == 0 {
            textraw[..iwidth].fill(b' ');
            textraw[iwidth] = b'\n';
        }
        if lenod > 0 {
            textraw[..12].copy_from_slice(b"C_SU_SETID,\0");
            textraw[11..11 + Rod.len()].copy_from_slice(Rod.as_bytes());
            textraw[11 + Rod.len()] = 0;
            if iwtype == 0 {
                textraw[11 + lenod] = b' ';
            } else {
                textraw[11 + lenod] = b'\n';
                textraw[12 + lenod] = 0;
            }
        } else if isetid == 1 {
            textraw[..15].copy_from_slice(b"C_SU_SETID,ANY\0");
            if iwtype == 0 {
                textraw[14] = b' ';
            } else {
                textraw[14] = b'\n';
                textraw[15] = 0;
            }
        } else {
            textraw[..16].copy_from_slice(b"C_SU_SETID,NONE\0");
            if iwtype == 0 {
                textraw[15] = b' ';
            } else {
                textraw[15] = b'\n';
                textraw[16] = 0;
            }
        }
        let end = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
        fputs_bytes(&textraw[..end], &fpW);

        let mut ibeg = if isetid == 0 { 0 } else { 1 };
        let mut iend = num_names;
        if iwtype != 0 {
            ibeg = 0;
            iend = numcases;
        }
        // Both C blocks write into the same initialized textraw array.
        // Filling only iwidth bytes preserves suffixes from longer records.
        for (title, fields) in [("C_SU_FORMS", &form2), ("C_SU_NAMES", &name2)] {
            if iwtype == 0 {
                textraw[..iwidth].fill(b' ');
                textraw[iwidth] = b'\n';
            }
            textraw[..10].copy_from_slice(title.as_bytes());
            textraw[10] = if iwtype == 0 { b' ' } else { b'\n' };
            if iwtype != 0 {
                textraw[11] = 0;
            }
            let end = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
            fputs_bytes(&textraw[..end], &fpW);

            if iwtype == 0 {
                textraw[..iwidth].fill(b' ');
                textraw[iwidth] = b'\n';
            }
            mspot = 0;
            if isetid == 1 {
                mspot = 8;
                textraw[..9].copy_from_slice(b"C_SU_ID,\0");
            }
            for i in ibeg..iend {
                let mleng = fields[i].len();
                if iwtype == 0 && mspot + mleng + 1 > iwidth {
                    textraw[mspot - 1] = b' ';
                    let end = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
                    fputs_bytes(&textraw[..end], &fpW);
                    textraw[..iwidth].fill(b' ');
                    textraw[iwidth] = b'\n';
                    mspot = 10;
                    textraw[..11].copy_from_slice(b"C_SU_MORE,\0");
                }
                textraw[mspot..mspot + mleng].copy_from_slice(fields[i].as_bytes());
                mspot += mleng;
                textraw[mspot] = b',';
                mspot += 1;
            }
            if iwtype == 0 {
                textraw[mspot - 1] = b' ';
            } else {
                textraw[mspot - 1] = b'\n';
                textraw[mspot] = 0;
            }
            let end = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
            fputs_bytes(&textraw[..end], &fpW);
        }
    }

    let mut count = 0usize;
    let mut mcount = 0i32;
    let mut ncount = 0i32;
    let mut lenerr = 0i32;
    let mut comerr = 0i32;
    let mut morerr = 0i32;
    let mut numerr = 0i32;
    let mut lrgerr = 0i32;
    let mut nblank = 0i32;
    let mut nextrow = 0;
    let mut nup = 0;
    let mut nlfield = -999999999999999i64;
    let incint = 100000000000000i64;
    let mut nrep = incint;

    loop {
        // Like fgets, retain bytes past the terminating NUL from the previous
        // record, and split long physical lines after maxtext - 1 bytes.
        let mut nread = 0;
        while nread < maxtext - 1 {
            let c = fgetc(&fpR);
            if c == EOF {
                break;
            }
            textraw[nread] = c as u8;
            nread += 1;
            if c == b'\n' as i32 {
                break;
            }
        }
        if nread == 0 {
            break;
        }
        textraw[nread] = 0;
        ncount += 1;
        if ncount >= nicerecord {
            let mut textfront = textraw[..10].to_vec();
            textfront.make_ascii_lowercase();
            if textfront.starts_with(b"c_su") || nextrow == 1 {
                nextrow = 0;
                if textfront.starts_with(b"c_su_names") || textfront.starts_with(b"c_su_forms") {
                    nextrow = 1;
                }
            } else if lenid < 1 || textraw.starts_with(&Rid.as_bytes()[..lenid]) {
                if count < numR {
                    let lenraw = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
                    if irtype == 0 {
                        for n in 0..pross.len() {
                            let lead = (plead[n] - 1) as usize;
                            let trail = ptrail[n] as usize;
                            let plen = trail - lead;
                            // strncpy pads with zero after the first NUL, and leaves
                            // the rest of textbeg from the preceding operation.
                            let mut ended = false;
                            for m in 0..plen {
                                let b = textraw[lead + m];
                                if b == 0 {
                                    ended = true;
                                }
                                textbeg[m] = if ended { 0 } else { b };
                            }
                            if pflag[n] == 2 || pflag[n] == 3 {
                                let mut mlead = 0;
                                for m in 0..plen {
                                    if (b'1'..=b'9').contains(&textbeg[m])
                                        || matches!(textbeg[m], b'+' | b'-')
                                    {
                                        mlead = m;
                                        if matches!(textbeg[m], b'+' | b'-') {
                                            mlead += 1;
                                        }
                                        break;
                                    } else {
                                        textbeg[m] = b' ';
                                    }
                                }
                                let mut mtrail = plen - 1;
                                for m in (0..plen).rev() {
                                    if textbeg[m].is_ascii_digit() {
                                        mtrail = m;
                                        break;
                                    } else {
                                        textbeg[m] = b' ';
                                    }
                                }
                                if pflag[n] == 3 {
                                    for m in mlead..mtrail {
                                        if !textbeg[m].is_ascii_digit() {
                                            textbeg[m] = b'0';
                                        }
                                    }
                                }
                                let mut ended = false;
                                for m in 0..plen {
                                    let b = textbeg[m];
                                    if b == 0 {
                                        ended = true;
                                    }
                                    textraw[lead + m] = if ended { 0 } else { b };
                                }
                            } else if pflag[n] == -1 || pflag[n] == -2 {
                                let end = textbeg.iter().position(|b| *b == 0).unwrap_or(maxtext);
                                let (mut dv, used) =
                                    strtod(&String::from_utf8_lossy(&textbeg[..end]));
                                if used > 0 {
                                    if pflag[n] == -1 {
                                        dv += pvalu[n];
                                    } else {
                                        dv *= pvalu[n];
                                    }
                                    let formatted = sprintf!("%.20f", dv);
                                    textraw2[..formatted.len()]
                                        .copy_from_slice(formatted.as_bytes());
                                    textraw2[formatted.len()] = 0;
                                    let mut j = 101;
                                    for m in 0..101 {
                                        if textraw2[m] != b' ' {
                                            j = m;
                                            break;
                                        }
                                    }
                                    let mut ended = false;
                                    for m in 0..plen {
                                        let b = textraw2[j + m];
                                        if b == 0 {
                                            ended = true;
                                        }
                                        textraw[lead + m] = if ended { 0 } else { b };
                                    }
                                } else {
                                    textraw[lead..trail].fill(b' ');
                                    textraw[trail - 1] = b'*';
                                }
                            } else if pflag[n] == 0 {
                                textraw[lead..trail].fill(b'0');
                            } else if pflag[n] == 1 {
                                textraw[lead..trail].fill(b' ');
                            }
                        }
                        let mut icomma = 0usize;
                        for ineed in 0..numcases {
                            if itrail[ineed] >= lenraw as i32 {
                                lenerr += 1;
                                if lenerr < 4 {
                                    warn!("Error at record %d   Record-too-short for requested fixed ranges", ncount);
                                    if lenerr == 3 {
                                        warn!("Have 3 Record-too-short warnings, no more will be printed.");
                                    }
                                }
                                break;
                            }
                            let lead = (ilead[ineed] - 1) as usize;
                            let plen = (itrail[ineed] - ilead[ineed] + 1) as usize;
                            let mut ended = false;
                            for m in 0..plen {
                                let b = textraw[lead + m];
                                if b == 0 {
                                    ended = true;
                                }
                                textbeg[icomma + m] = if ended { 0 } else { b };
                            }
                            icomma += plen;
                            textbeg[icomma] = rdel;
                            icomma += 1;
                        }
                        let mut ended = false;
                        for m in 0..icomma {
                            let b = textbeg[m];
                            if b == 0 {
                                ended = true;
                            }
                            textraw[m] = if ended { 0 } else { b };
                        }
                        textraw[icomma] = 0;
                    }
                    if iwchop > 0 {
                        let last = textraw.iter().position(|b| *b == 0).unwrap_or(maxtext);
                        if lenod > 0 {
                            textraw2[..Rod.len()].copy_from_slice(Rod.as_bytes());
                            textraw2[Rod.len()] = 0;
                            textraw2[lenod] = b',';
                            textraw2[lenod + 1..last + lenod + 1].copy_from_slice(&textraw[..last]);
                            textraw2[last + lenod] = b'\n';
                            textraw2[last + lenod + 1] = 0;
                        } else {
                            textraw2[..last].copy_from_slice(&textraw[..last]);
                            if last > 0 {
                                textraw2[last - 1] = b'\n';
                            }
                            textraw2[last] = 0;
                        }
                        let end = textraw2.iter().position(|b| *b == 0).unwrap_or(maxtext);
                        fputs_bytes(&textraw2[..end], &fpW);
                        if maxrecords == -2 {
                            mcount += 1;
                            continue;
                        }
                    }
                    getCSV(
                        &textraw,
                        &mut textbeg,
                        rdel,
                        &mut RecInfo[count].dfield,
                        &nspot[..numcases],
                        ncount,
                        &mut comerr,
                        &mut morerr,
                        &mut numerr,
                        &mut nblank,
                    );
                    for j in 0..numcases {
                        let v = RecInfo[count].dfield[j];
                        if v < 1.0e308 && (v > valmx[j] || 0. - v > valmx[j]) {
                            lrgerr += 1;
                            if lrgerr < 4 {
                                warn!(
                                    "Error at record %d number-too-large for SU name (%.2f)",
                                    ncount, v
                                );
                                if lrgerr == 3 {
                                    warn!("Have 3 number-too-large-for-SU-name warnings, no more will be printed.");
                                }
                            }
                        }
                    }
                    if mapx[3] > -1 {
                        let d = &mut RecInfo[count].dfield;
                        let cfr = d[mapx[0] as usize];
                        let ctr = d[mapx[1] as usize];
                        let cir = d[mapx[2] as usize];
                        let rfr = d[mapx[3] as usize];
                        let rtr = d[mapx[4] as usize];
                        if cfr > ctr {
                            d[mapx[0] as usize] = ctr;
                            d[mapx[1] as usize] = cfr;
                            d[mapx[3] as usize] = rtr;
                            d[mapx[4] as usize] = rfr;
                        }
                        if (d[mapx[0] as usize] - d[mapx[1] as usize]).abs() < dtolh * 2. {
                            d[mapx[2] as usize] = 1.;
                        } else {
                            d[mapx[2] as usize] = cir.abs();
                        }
                    }
                    for k in 0..num_to_sort_by {
                        RecInfo[count].lfield[k] =
                            longt(RecInfo[count].dfield[ktol[k]], dtolh, dtol);
                    }
                    if unrepeat > -2147483645 {
                        if nup > 0 {
                            if nlfield > RecInfo[count].lfield[0] {
                                nup = -1;
                                nrep = nrep.wrapping_add(incint);
                            }
                        } else if nup < 0 {
                            if nlfield < RecInfo[count].lfield[0] {
                                nup = 1;
                                nrep = nrep.wrapping_add(incint);
                            }
                        } else if nlfield > -999999999999999i64 {
                            nup = 1;
                            if nlfield > RecInfo[count].lfield[0] {
                                nup = -1;
                            }
                        }
                        nlfield = RecInfo[count].lfield[0];
                        RecInfo[count].lfield[0] = RecInfo[count].lfield[0].wrapping_add(nrep);
                    }
                    if iwchop > 0 {
                        if maxrecords > -1 {
                            count += 1;
                        }
                        mcount += 1;
                        continue;
                    }
                    if iwtype == 0 {
                        textraw2[..iwidth].fill(b' ');
                        textraw2[iwidth] = b'\n';
                        textraw2[..lenod].copy_from_slice(&Rod.as_bytes()[..lenod]);
                        for ineed in 0..numcases {
                            // C sutoolcsv.c:2652/2655 writes before textraw2
                            // when bare names leave column bounds at -1.
                            // Match the bounds check in the UB-fixed C.
                            if ilead[ineed] < 1
                                || itrail[ineed] < ilead[ineed]
                                || itrail[ineed] > iwidth as i32
                            {
                                err!("**** Error: fixed output requires column bounds in names=.");
                            }
                            let formatted = sprintf!(&forms[ineed], RecInfo[count].dfield[ineed]);
                            let mfill = formatted.len();
                            // C sprintf writes into textbeg. Later process=
                            // operations reuse its suffix after a short strncpy.
                            textbeg[..mfill].copy_from_slice(formatted.as_bytes());
                            textbeg[mfill] = 0;
                            if itrail[ineed] - mfill as i32 + 1 < ilead[ineed] {
                                textraw2[(itrail[ineed] - 1) as usize] = b'*';
                            } else {
                                let start = (itrail[ineed] - mfill as i32) as usize;
                                textraw2[start..start + mfill]
                                    .copy_from_slice(formatted.as_bytes());
                            }
                        }
                        let end = textraw2.iter().position(|b| *b == 0).unwrap_or(maxtext);
                        fputs_bytes(&textraw2[..end], &fpW);
                    } else {
                        textraw2[..lenod].copy_from_slice(&Rod.as_bytes()[..lenod]);
                        let mut mhere = lenod;
                        for ineed in 0..numcases {
                            let v = RecInfo[count].dfield[ineed];
                            if v < 1.0e308 {
                                let formatted = sprintf!(&forms[ineed], v);
                                let mfill = formatted.len();
                                textbeg[..mfill].copy_from_slice(formatted.as_bytes());
                                textbeg[mfill] = 0;
                                textraw2[mhere] = b',';
                                textraw2[mhere + 1..mhere + 1 + mfill]
                                    .copy_from_slice(formatted.as_bytes());
                                mhere += mfill + 1;
                            } else {
                                textraw2[mhere] = b',';
                                textraw2[mhere + 1] = b'*';
                                mhere += 2;
                            }
                        }
                        if isetid == 0 {
                            textraw2[0] = b' ';
                        }
                        textraw2[mhere] = b'\n';
                        textraw2[mhere + 1] = 0;
                        fputs_bytes(&textraw2[..mhere + 1], &fpW);
                    }
                }
                if maxrecords > -1 {
                    count += 1;
                }
                mcount += 1;
            }
        }
    }
    if nblank > 0 {
        warn!("Total all-blank fields: %d. Assumed zero for all.", nblank);
    }
    if numerr > 0 {
        warn!(
            "Total Field-unreadable as a number:        %d (will error-halt SUGEOMCSV)",
            numerr
        );
    }
    if morerr > 0 {
        warn!(
            "Total Two-numbers in one field:            %d (will error-halt SUGEOMCSV)",
            morerr
        );
    }
    if lenerr > 0 {
        warn!(
            "Total Record-too-short to get all values:  %d (will error-halt SUGEOMCSV)",
            lenerr
        );
    }
    if comerr > 0 {
        warn!(
            "Total Not-enough-commas to get all values: %d (will error-halt SUGEOMCSV)",
            comerr
        );
    }
    if lrgerr > 0 {
        warn!(
            "Total Number-too-large-for-SU-name:        %d (will error-halt SUGEOMCSV)",
            lrgerr
        );
    }
    if unrepeat > -2147483645 {
        warn!(
            "For unrepeat option, the text incrementing integer ended at: %d ",
            nrep / incint
        );
    }
    warn!(
        "Have allocated memory to store values from %d records. Found %d records.",
        numR as i32, mcount
    );
    if maxrecords > -1 && mcount > numR as i32 {
        warn!(
            "Error: Too many records read (%d) for your maxrecords= value (%d).",
            mcount, numR as i32
        );
        // The vendor proceeds into qsort with more entries than allocated
        // and dies with SIGSEGV before its buffered output file is flushed.
        vendor_segv();
    }
    if maxrecords < 0 {
        exit(0);
    }
    let numR = mcount as usize;
    RecInfo.truncate(numR);
    // C uses qsort(RecInfo, numR, sizeof(struct PointInfo), compSort).
    let mut RecInfo: Vec<PointView<'_>> = RecInfo
        .iter()
        .map(|point| PointView {
            dfield: &point.dfield,
            lfield: &point.lfield,
        })
        .collect();
    // SAFETY: qsort receives the initialized, exclusively borrowed view array
    // and its exact element size. Each view is Copy and contains references
    // that remain valid while qsort rearranges the views. compSort only reads
    // these views and their equally sized key arrays; it cannot unwind.
    unsafe {
        qsort(
            RecInfo.as_mut_ptr().cast(),
            numR,
            std::mem::size_of::<PointView<'_>>(),
            compSort,
        );
    }
    checks::check(&RecInfo, &mapx, &match_, dtolh, dtol);
    if numerr > 0 || morerr > 0 || lenerr > 0 || comerr > 0 {
        warn!("File has Field-unreadable, Two-numbers, or Record-short errors. These often cause multiple subsequent error/warnings.");
    }
    exit(0);
}

// The C file's private helpers.
unsafe extern "C" fn compSort(q1: *const std::ffi::c_void, q2: *const std::ffi::c_void) -> i32 {
    // SAFETY: qsort passes pointers to initialized PointView elements (or
    // temporary copies of them), aligned to the element's pointer alignment.
    let p1 = unsafe { &*q1.cast::<PointView<'_>>() };
    let p2 = unsafe { &*q2.cast::<PointView<'_>>() };
    for n in 0..p1.lfield.len() {
        if p1.lfield[n] < p2.lfield[n] {
            return -1;
        }
        if p1.lfield[n] > p2.lfield[n] {
            return 1;
        }
    }
    0
}

// Preserve the vendor's observed SIGSEGV status on malformed inputs without
// invalid memory accesses or direct C calls. exec replaces this process, so
// the shell's default signal disposition also avoids Rust's signal handler.
// No buffered output is flushed, just as on the vendor's fatal signal.
fn vendor_segv() -> ! {
    use std::os::unix::process::CommandExt;
    let _ = std::process::Command::new("/bin/sh")
        .args(["-c", "kill -s SEGV $$"])
        .exec();
    exit(139);
}

fn countRec(rfile: &CFile, textraw: &mut [u8], Rid: &str, lenid: usize, nicerecord: i32) -> usize {
    let mut count = 0;
    let mut ncount = 0;
    let mut nextrow = 0;
    loop {
        let mut nread = 0;
        while nread < textraw.len() - 1 {
            let c = fgetc(rfile);
            if c == EOF {
                break;
            }
            textraw[nread] = c as u8;
            nread += 1;
            if c == b'\n' as i32 {
                break;
            }
        }
        if nread == 0 {
            break;
        }
        textraw[nread] = 0;
        ncount += 1;
        if ncount >= nicerecord {
            let mut textfront = textraw[..10].to_vec();
            textfront.make_ascii_lowercase();
            if textfront.starts_with(b"c_su") || nextrow == 1 {
                nextrow = 0;
                if textfront.starts_with(b"c_su_names") || textfront.starts_with(b"c_su_forms") {
                    nextrow = 1;
                }
            } else if lenid < 1 || textraw.starts_with(&Rid.as_bytes()[..lenid]) {
                count += 1;
            }
        }
    }
    fseeko(rfile, 0, SEEK_SET);
    count
}

#[allow(clippy::too_many_arguments)]
fn getCSV(
    textraw: &[u8],
    textbeg: &mut [u8],
    rdel: u8,
    dfield: &mut [f64],
    nspot: &[usize],
    ncount: i32,
    comerr: &mut i32,
    morerr: &mut i32,
    numerr: &mut i32,
    nblank: &mut i32,
) {
    let mut nbeg = -1i32;
    let mut nfield = 0;
    let mut ineed = 0;
    for n in 0..textraw.len() {
        if textraw[n] == rdel || matches!(textraw[n], 0 | b'\n' | b'\r') {
            if nfield == nspot[ineed] {
                let mut dval = 1.1e308;
                let mut nb = -1;
                let plen = (n as i32 - nbeg - 1) as usize;
                if plen > 0 {
                    textbeg[..plen].copy_from_slice(&textraw[(nbeg + 1) as usize..n]);
                    textbeg[plen] = 0;
                    let mut ib = -1;
                    for m in 0..plen {
                        if textbeg[m] != b' ' {
                            nb = m as i32;
                            if ib > -1 {
                                *morerr += 1;
                                if *morerr < 4 {
                                    let message = sprintf!(
                                        "Error at record %d   two-numbers in field (%s)",
                                        ncount,
                                        ""
                                    );
                                    warn_field(&message, &textbeg[..plen]);
                                    if *morerr == 3 {
                                        warn!("Have 3 two-numbers in field warnings, no more will be printed.");
                                    }
                                }
                                nb = -1;
                                break;
                            }
                        }
                        if textbeg[m] == b' ' && nb > -1 {
                            ib = m as i32;
                        }
                    }
                    if nb > -1 {
                        let end = textbeg[..plen].iter().position(|b| *b == 0).unwrap_or(plen);
                        let (value, used) = strtod(&String::from_utf8_lossy(&textbeg[..end]));
                        if used > 0 {
                            dval = value;
                        } else {
                            *numerr += 1;
                            if *numerr < 4 {
                                let message = sprintf!(
                                    "Error at record %d   field-unreadable as a number (%s)",
                                    ncount,
                                    ""
                                );
                                warn_field(&message, &textbeg[..end]);
                                if *numerr == 3 {
                                    warn!("Have 3 field-unreadable warnings, no more will be printed.");
                                }
                            }
                        }
                    }
                }
                if nb < 0 {
                    *nblank += 1;
                    dval = 0.;
                }
                dfield[ineed] = dval;
                ineed += 1;
                if ineed >= nspot.len() {
                    break;
                }
            }
            if matches!(textraw[n], 0 | b'\n' | b'\r') && ineed < nspot.len() {
                *comerr += 1;
                if *comerr < 4 {
                    warn!(
                        "Error at record %d   Not-enough-commas in record to get all values",
                        ncount
                    );
                    if *comerr == 3 {
                        warn!("Have 3 Not-enough-comma warnings, no more will be printed.");
                    }
                }
            }
            nbeg = n as i32;
            nfield += 1;
        }
    }
}

// The two getCSV field diagnostics need a raw byte insertion before their
// closing ')'. Core warn!/sprintf! accept UTF-8 strings only; use core's
// formatting and byte output primitives without decoding the field.
fn warn_field(message: &str, field: &[u8]) {
    let end = field.iter().position(|b| *b == 0).unwrap_or(field.len());
    if let Ok(field) = std::str::from_utf8(&field[..end]) {
        warn!("%s%s)", &message[..message.len() - 1], field);
    } else {
        if fflush(&stdout()) == EOF {
            fputs("\nwarn: fflush failed on stdout", &stderr());
        }
        fprintf!(
            stderr(),
            "\n%s: %s",
            xargv().first().map(String::as_str).unwrap_or(""),
            &message[..message.len() - 1]
        );
        fputs_bytes(&field[..end], &stderr());
        fputs(")\n", &stderr());
    }
}

pub(crate) fn longt(dvalue: f64, dtolh: f64, dtol: f64) -> i64 {
    if dvalue >= 0.0 {
        ((dvalue + dtolh) * dtol) as i64
    } else {
        ((dvalue - dtolh) * dtol) as i64
    }
}

fn tparse(tbuf: &[u8], d: u8) -> Vec<String> {
    tbuf.split(|b| *b == d)
        .map(|s| {
            if s.is_empty() {
                "null".into()
            } else {
                String::from_utf8_lossy(s).into_owned()
            }
        })
        .collect()
}

fn GetCase(cbuf: &str) -> i32 {
    if cbuf.starts_with("null") {
        return 0;
    }
    match cbuf {
        "tracl" => 1,
        "tracr" => 2,
        "fldr" => 3,
        "tracf" => 4,
        "ep" => 5,
        "cdp" => 6,
        "cdpt" => 7,
        "trid" => 8,
        "nvs" => 9,
        "nhs" => 10,
        "duse" => 11,
        "offset" => 12,
        "gelev" => 13,
        "selev" => 14,
        "sdepth" => 15,
        "gdel" => 16,
        "sdel" => 17,
        "swdep" => 18,
        "gwdep" => 19,
        "scalel" => 20,
        "scalco" => 21,
        "sx" => 22,
        "sy" => 23,
        "gx" => 24,
        "gy" => 25,
        "counit" => 26,
        "wevel" => 27,
        "swevel" => 28,
        "sut" => 29,
        "gut" => 30,
        "sstat" => 31,
        "gstat" => 32,
        "tstat" => 33,
        "laga" => 34,
        "lagb" => 35,
        "delrt" => 36,
        "muts" => 37,
        "mute" => 38,
        "ns" => 39,
        "dt" => 40,
        "gain" => 41,
        "igc" => 42,
        "igi" => 43,
        "corr" => 44,
        "sfs" => 45,
        "sfe" => 46,
        "slen" => 47,
        "styp" => 48,
        "stas" => 49,
        "stae" => 50,
        "tatyp" => 51,
        "afilf" => 52,
        "afils" => 53,
        "nofilf" => 54,
        "nofils" => 55,
        "lcf" => 56,
        "hcf" => 57,
        "lcs" => 58,
        "hcs" => 59,
        "year" => 60,
        "day" => 61,
        "hour" => 62,
        "minute" => 63,
        "sec" => 64,
        "timbas" => 65,
        "trwf" => 66,
        "grnors" => 67,
        "grnofr" => 68,
        "grnlof" => 69,
        "gaps" => 70,
        "otrav" => 71,
        "d1" => 72,
        "f1" => 73,
        "d2" => 74,
        "f2" => 75,
        "ungpow" => 76,
        "unscale" => 77,
        "ntr" => 78,
        "mark" => 79,
        _ if cbuf.starts_with("numb") => 1000i32.wrapping_add(atoi(&cbuf[4..])),
        _ => -1,
    }
}
