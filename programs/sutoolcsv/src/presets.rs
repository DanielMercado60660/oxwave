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

//! Literal SPS name and format expansions from su/main/headers/sutoolcsv.c.
//! Credits: Andre Latour; started from sugeom by Fernando Roxo <petro@roxo.org>.

pub fn names(version: i32, id: &str, all: bool) -> Option<Vec<String>> {
    let a: &[&str] = match (version, id, all) {
        (2, "X", false) => &[
            "c_su_id",
            "2_null2_7",
            "8_match1_15",
            "16_null4_16",
            "17_null5_17",
            "18_grnofr_27",
            "28_grnlof_37",
            "38_null8_38",
            "39_matche1_cf_43",
            "44_matche1_ct_48",
            "49_matche1_ci_49",
            "50_grnors_59",
            "60_gaps_rf_69",
            "70_gaps_rt_79",
            "80_null15_80",
        ],
        (2, "X", true) => &[
            "c_su_id",
            "2_numb2_7",
            "8_match1_15",
            "16_numb4_16",
            "17_numb5_17",
            "18_grnofr_27",
            "28_grnlof_37",
            "38_numb8_38",
            "39_matche1_cf_43",
            "44_matche1_ct_48",
            "49_matche1_ci_49",
            "50_grnors_59",
            "60_gaps_rf_69",
            "70_gaps_rt_79",
            "80_numb15_80",
        ],
        (2, "S", false) => &[
            "c_su_id",
            "2_grnofr_11",
            "12_grnlof_21",
            "22_null4_23",
            "24_null5_24",
            "25_null6_26",
            "27_sstat_30",
            "31_sdepth_34",
            "35_sdel_38",
            "39_sut_40",
            "41_swdep_46",
            "47_sx_55",
            "56_sy_65",
            "66_selev_71",
            "72_null15_74",
            "75_null16_76",
            "77_null17_78",
            "79_null18_80",
        ],
        (2, "S", true) => &[
            "c_su_id",
            "2_grnofr_11",
            "12_grnlof_21",
            "22_numb4_23",
            "24_numb5_24",
            "25_numb6_26",
            "27_sstat_30",
            "31_sdepth_34",
            "35_sdel_38",
            "39_sut_40",
            "41_swdep_46",
            "47_sx_55",
            "56_sy_65",
            "66_selev_71",
            "72_numb15_74",
            "75_numb16_76",
            "77_numb17_78",
            "79_numb18_80",
        ],
        (2, "R", false) => &[
            "c_su_id",
            "2_grnors_11",
            "12_gaps_21",
            "22_null4_23",
            "24_null5_24",
            "25_null6_26",
            "27_gstat_30",
            "31_null8_34",
            "35_gdel_38",
            "39_gut_40",
            "41_gwdep_46",
            "47_gx_55",
            "56_gy_65",
            "66_gelev_71",
            "72_null15_74",
            "75_null16_76",
            "77_null17_78",
            "79_null18_80",
        ],
        (2, "R", true) => &[
            "c_su_id",
            "2_grnors_11",
            "12_gaps_21",
            "22_numb4_23",
            "24_numb5_24",
            "25_numb6_26",
            "27_gstat_30",
            "31_numb8_34",
            "35_gdel_38",
            "39_gut_40",
            "41_gwdep_46",
            "47_gx_55",
            "56_gy_65",
            "66_gelev_71",
            "72_numb15_74",
            "75_numb16_76",
            "77_numb17_78",
            "79_numb18_80",
        ],
        (1, "X", false) => &[
            "c_su_id",
            "2_null2_7",
            "8_match1_11",
            "12_null4_12",
            "13_null5_13",
            "14_grnofr_29",
            "30_grnlof_37",
            "38_null8_38",
            "39_matche1_cf_42",
            "43_matche1_ct_46",
            "47_matche1_ci_47",
            "48_grnors_63",
            "64_gaps_rf_71",
            "72_gaps_rt_79",
            "80_null15_80",
        ],
        (1, "X", true) => &[
            "c_su_id",
            "2_numb2_7",
            "8_match1_11",
            "12_numb4_12",
            "13_numb5_13",
            "14_grnofr_29",
            "30_grnlof_37",
            "38_numb8_38",
            "39_matche1_cf_42",
            "43_matche1_ct_46",
            "47_matche1_ci_47",
            "48_grnors_63",
            "64_gaps_rf_71",
            "72_gaps_rt_79",
            "80_numb15_80",
        ],
        (1, "S", false) => &[
            "c_su_id",
            "2_grnofr_17",
            "18_grnlof_25",
            "26_null4_26",
            "26_null5_26",
            "27_null6_28",
            "29_sstat_32",
            "33_sdepth_36",
            "37_sdel_40",
            "41_sut_42",
            "43_swdep_46",
            "47_sx_55",
            "56_sy_65",
            "66_selev_71",
            "72_null15_74",
            "75_null16_76",
            "77_null17_78",
            "79_null18_80",
        ],
        (1, "S", true) => &[
            "c_su_id",
            "2_grnofr_17",
            "18_grnlof_25",
            "26_numb4_26",
            "26_numb5_26",
            "27_numb6_28",
            "29_sstat_32",
            "33_sdepth_36",
            "37_sdel_40",
            "41_sut_42",
            "43_swdep_46",
            "47_sx_55",
            "56_sy_65",
            "66_selev_71",
            "72_numb15_74",
            "75_numb16_76",
            "77_numb17_78",
            "79_numb18_80",
        ],
        (1, "R", false) => &[
            "c_su_id",
            "2_grnors_17",
            "18_gaps_25",
            "26_null4_26",
            "26_null5_26",
            "27_null6_28",
            "29_gstat_32",
            "33_null8_36",
            "37_gdel_40",
            "41_gut_42",
            "43_gwdep_46",
            "47_gx_55",
            "56_gy_65",
            "66_gelev_71",
            "72_null15_74",
            "75_null16_76",
            "77_null17_78",
            "79_null18_80",
        ],
        (1, "R", true) => &[
            "c_su_id",
            "2_grnors_17",
            "18_gaps_25",
            "26_numb4_26",
            "26_numb5_26",
            "27_numb6_28",
            "29_gstat_32",
            "33_numb8_36",
            "37_gdel_40",
            "41_gut_42",
            "43_gwdep_46",
            "47_gx_55",
            "56_gy_65",
            "66_gelev_71",
            "72_numb15_74",
            "75_numb16_76",
            "77_numb17_78",
            "79_numb18_80",
        ],
        _ => return None,
    };
    Some(a.iter().map(|s| s.to_string()).collect())
}

pub fn forms(id: &str) -> Option<Vec<String>> {
    let a: &[&str] = match id {
        "X" => &[
            "%s", "%.0f", "%.0f", "%.0f", "%.0f", "%.2f", "%.2f", "%.0f", "%.0f", "%.0f", "%.0f",
            "%.2f", "%.2f", "%.2f", "%.0f",
        ],
        "S" => &[
            "%s", "%.2f", "%.2f", "%.0f", "%.0f", "%.0f", "%.0f", "%.1f", "%.0f", "%.0f", "%.1f",
            "%.1f", "%.1f", "%.1f", "%.0f", "%.0f", "%.0f", "%.0f",
        ],
        "R" => &[
            "%s", "%.2f", "%.2f", "%.0f", "%.0f", "%.0f", "%.0f", "%.1f", "%.0f", "%.0f", "%.1f",
            "%.1f", "%.1f", "%.1f", "%.0f", "%.0f", "%.0f", "%.0f",
        ],
        _ => return None,
    };
    Some(a.iter().map(|s| s.to_string()).collect())
}
