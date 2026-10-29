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

//! Port of comp/dct/lib/getfilter.c, wavepack1.c, wavepack2.c,
//! wavetrans1.c and wavetrans2.c.
//! Credits: Tong Chen, 05/18/94, 05/25/94, 04/05/95 and 04/12/95.

/// Owned C waveFilter: low-pass and quadrature mirror high-pass coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct WaveFilter {
    pub id: i32,
    pub len: usize,
    pub ishift: i32,
    pub filterh: Vec<f32>,
    pub filterg: Vec<f32>,
}

const WAVE_TAB: [&[f32]; 15] = [
    &[0.707106781, 0.707106781],
    &[
        0.4829629131445341,
        0.8365163037378077,
        0.2241438680420134,
        -0.1294095225512603,
    ],
    &[
        0.3326705529500825,
        0.8068915093110924,
        0.4598775021184914,
        -0.1350110200102546,
        -0.0854412738820267,
        0.0352262918857095,
    ],
    &[
        0.2303778133088964,
        0.7148465705529154,
        0.6308807679298587,
        -0.0279837694168599,
        -0.1870348117190931,
        0.0308413818355607,
        0.0328830116668852,
        -0.0105974017850690,
    ],
    &[
        0.1601023979741929,
        0.6038292697971895,
        0.7243085284377726,
        0.1384281459013203,
        -0.2422948870663823,
        -0.0322448695846381,
        0.0775714938400459,
        -0.0062414902127983,
        -0.0125807519990820,
        0.0033357252854738,
    ],
    &[
        0.1115407433501095,
        0.4946238903984533,
        0.7511339080210959,
        0.3152503517091982,
        -0.2262646939654400,
        -0.1297668675672625,
        0.0975016055873225,
        0.0275228655303053,
        -0.0315820393174862,
        0.0005538422011614,
        0.0047772575109455,
        -0.0010773010853085,
    ],
    &[
        0.0778520540850037,
        0.3965393194818912,
        0.7291320908461957,
        0.4697822874051889,
        -0.1439060039285212,
        -0.2240361849938412,
        0.0713092192668272,
        0.0806126091510774,
        -0.0380299369350104,
        -0.0165745416306655,
        0.0125509985560986,
        0.0004295779729214,
        -0.0018016407040473,
        0.0003537137999745,
    ],
    &[
        0.0544158422431072,
        0.3128715909143166,
        0.6756307362973195,
        0.5853546836542159,
        -0.0158291052563823,
        -0.2840155429615824,
        0.0004724845739124,
        0.1287474266204893,
        -0.017369301001809,
        -0.0440882539307971,
        0.0139810279174001,
        0.0087460940474065,
        -0.0048703529934520,
        -0.0003917403733770,
        0.0006754494064506,
        -0.0001174767841248,
    ],
    &[
        0.0380779473638778,
        0.2438346746125858,
        0.6048231236900955,
        0.6572880780512736,
        0.1331973858249883,
        -0.2932737832791663,
        -0.0968407832229492,
        0.1485407493381256,
        0.0307256814793385,
        -0.0676328290613279,
        0.0002509471148340,
        0.0223616621236798,
        -0.0047232047577518,
        -0.0042815036824635,
        0.0018476468830563,
        0.0002303857635232,
        -0.0002519631889427,
        0.0000393473203163,
    ],
    &[
        0.0266700579005473,
        0.1881768000776347,
        0.5272011889315757,
        0.6884590394534363,
        0.2811723436605715,
        -0.2498464243271598,
        -0.1959462743772862,
        0.1273693403357541,
        0.0930573646035547,
        -0.0713941471663501,
        -0.0294575368218399,
        0.0332126740593612,
        0.0036065535669870,
        -0.0107331754833007,
        0.0013953517470688,
        0.0019924052951925,
        -0.0006858566949564,
        -0.0001164668551285,
        0.0000935886703202,
        -0.0000132642028945,
    ],
    &[
        0.038580777747887,
        -0.12696912539621,
        -0.077161555495774,
        0.60749164138568,
        0.74568755893443,
        0.22658426519707,
    ],
    &[
        0.016387336464,
        -0.041464936782,
        -0.067372554722,
        0.386110066823,
        0.812723635450,
        0.416998396426,
        -0.076488599079,
        -0.059434418647,
        0.023680171946,
        0.005611434819,
        -0.001823208871,
        -0.000720549445,
    ],
    &[
        -0.0037935128643778759,
        0.0077825964256797869,
        0.023452696142119103,
        -0.065771911281431228,
        -0.061123390002955698,
        0.40517690240961679,
        0.79377722262562034,
        0.42848347637761869,
        -0.071799821619170590,
        -0.082301927106320283,
        0.034555027573344464,
        0.015880544863615901,
        -0.0090079761367322896,
        -0.0025745176881279692,
        0.0011175187708269618,
        0.00046621695982014403,
        -0.000070983302505704928,
        -0.000034599773197402695,
    ],
    &[
        0.00089231366822027571,
        -0.001629492013110849,
        -0.0073461663276562349,
        0.016068943964069236,
        0.026682300155628804,
        -0.081266699680313054,
        -0.05607731331647195,
        0.41530840703043015,
        0.78223893092049879,
        0.43438605649146839,
        -0.066627474263000752,
        -0.096220442033563697,
        0.039334427122913219,
        0.025082261845146933,
        -0.015211731527239149,
        -0.005658286685946038,
        0.0037514361569249027,
        0.0012665619286795187,
        -0.00058902075681143784,
        -0.00025997455231942175,
        0.000062339033865764618,
        0.000031229876078043358,
        -0.0000032596804448576129,
        0.0000017849845586999338,
    ],
    &[
        -0.00021208086333630681,
        0.0003585896772556986,
        0.0021782363048412847,
        -0.0041593587816039935,
        -0.010131117538045594,
        0.023408156761592795,
        0.028168029062141497,
        -0.091920010548806413,
        -0.05204316321623739,
        0.42156620672876544,
        0.77428960374028455,
        0.43799162622836413,
        -0.062035963905608969,
        -0.10557420870583534,
        0.041289208740734169,
        0.032683574283249535,
        -0.019761779011723959,
        -0.0091642311530462268,
        0.00676418541866332,
        0.0024333732092240538,
        -0.0016628637690858134,
        -0.00063813129615137752,
        0.00030225951979184068,
        0.00014054114890107723,
        -0.000041340484491956856,
        -0.000021315014062244917,
        0.0000037345967496715605,
        0.0000020638063902331633,
        -0.00000016740829374930063,
        -0.000000095157917046829356,
    ],
];

/// C prototype: `waveFilter *waveGetfilter(int id);`
/// Returns None for an unsupported id. Negative ids are undefined in C.
pub fn waveGetfilter(id: i32) -> Option<WaveFilter> {
    let h = WAVE_TAB.get(id as usize)?.to_vec();
    let len = h.len();
    let mut g = vec![0.0; len];
    for i in 0..len {
        g[i] = h[len - 1 - i];
    }
    for i in (1..len).step_by(2) {
        g[i] = -g[i];
    }
    let ishift = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -4, -6, -8, -10][id as usize];
    Some(WaveFilter {
        id,
        len,
        ishift,
        filterh: h,
        filterg: g,
    })
}

/// C prototype: `void waveFreefilter(waveFilter *filter);`
pub fn waveFreefilter(filter: WaveFilter) {
    drop(filter);
}

// The C row kernels repeat the 1D kernels with an innermost column loop.
// Keep that loop here; the 1D entry points use a single column. All f32
// additions and multiplications retain the C ordering.
fn position(n: usize, ishift: i32) -> usize {
    let mut pos = n;
    while pos < ishift.unsigned_abs() as usize {
        pos <<= 1;
    }
    (pos as i32 + ishift) as usize
}

fn pack_rows(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    filter: &WaveFilter,
    pow: usize,
    stage: usize,
    typ: i32,
) -> i32 {
    if stage > pow {
        return 0;
    }
    let n = 1 << pow;
    let nc = x[0].len();
    if stage == 0 {
        for i in 0..n {
            for l in 0..nc {
                if typ == 0 {
                    y[i][l] = x[i][l];
                } else {
                    x[i][l] = y[i][l];
                }
            }
        }
        return 1;
    }
    let len = filter.len;
    let h = &filter.filterh;
    let g = &filter.filterg;
    let nhalf = n >> 1;
    let mut tmp = vec![vec![0.0f32; nc]; n];
    if typ == 0 {
        let pos = position(n, filter.ishift);
        for j in 0..nhalf {
            for l in 0..nc {
                tmp[j][l] = 0.0;
                y[j][l] = 0.0;
            }
            for i in (1..len).step_by(2) {
                let k = (i + 2 * j - 1 + pos) % n;
                for l in 0..nc {
                    tmp[j][l] += x[k][l] * h[i - 1] + x[k + 1][l] * h[i];
                    y[j][l] += x[k][l] * g[i - 1] + x[k + 1][l] * g[i];
                }
            }
        }
        for ilevel in 1..stage {
            let nband = (1 << ilevel) >> 1;
            let nmid = 1 << (pow - ilevel);
            let nmidh = nmid >> 1;
            let pos = position(nmid, filter.ishift);
            for iband in 0..nband {
                let ifirst = (iband * nmid) << 1;
                let mut ofirst = ifirst + nmid;
                for j in 0..nmidh {
                    let jj = ofirst + j;
                    for l in 0..nc {
                        tmp[jj][l] = 0.0;
                        y[jj][l] = 0.0;
                    }
                    for i in (1..len).step_by(2) {
                        let k = (i + 2 * j - 1 + pos) % nmid + ifirst;
                        for l in 0..nc {
                            tmp[jj][l] += tmp[k][l] * h[i - 1] + tmp[k + 1][l] * h[i];
                            y[jj][l] += tmp[k][l] * g[i - 1] + tmp[k + 1][l] * g[i];
                        }
                    }
                }
                ofirst += nmidh;
                for j in 0..nmidh {
                    let jj = ofirst + j;
                    for l in 0..nc {
                        tmp[jj][l] = 0.0;
                        y[jj][l] = 0.0;
                    }
                    for i in (1..len).step_by(2) {
                        let k = (i + 2 * j - 1 + pos) % nmid + ifirst;
                        for l in 0..nc {
                            tmp[jj][l] += y[k][l] * g[i - 1] + y[k + 1][l] * g[i];
                            y[jj][l] += y[k][l] * h[i - 1] + y[k + 1][l] * h[i];
                        }
                    }
                }
                for j in 0..nmidh {
                    for l in 0..nc {
                        tmp[ifirst + j][l] = tmp[ofirst + j][l];
                        y[ifirst + j][l] = y[ofirst + j][l];
                    }
                }
            }
        }
        let nmid = 1 << (pow - stage);
        let nmidh = nmid << 1;
        let nband = (1 << stage) >> 1;
        for iband in 0..nband {
            let ifirst = iband * nmidh;
            let ofirst = ifirst + nmid;
            for i in 0..nmid {
                for l in 0..nc {
                    y[ofirst + i][l] = tmp[ifirst + i][l];
                }
            }
        }
    } else {
        for i in 0..n {
            for l in 0..nc {
                x[i][l] = y[i][l];
            }
        }
        for ilevel in (2..=stage).rev() {
            let nmid = 1 << (pow - ilevel);
            let nmidh = nmid << 1;
            let nband = (1 << ilevel) >> 1;
            let pos = position(nmidh, filter.ishift);
            for iband in 0..nband {
                let ifirst = iband * nmidh;
                let ofirst = ifirst + nmid;
                if iband == ((iband >> 1) << 1) {
                    for i in 0..nmid {
                        for l in 0..nc {
                            tmp[ifirst + i][l] = x[ofirst + i][l];
                        }
                    }
                    for i in 0..nmid {
                        for l in 0..nc {
                            tmp[ofirst + i][l] = x[ifirst + i][l];
                        }
                    }
                } else {
                    for i in 0..nmidh {
                        for l in 0..nc {
                            tmp[ifirst + i][l] = x[ifirst + i][l];
                        }
                    }
                }
                for i in 0..nmidh {
                    for l in 0..nc {
                        x[ifirst + i][l] = 0.0;
                    }
                }
                for j in 0..nmid {
                    let jj = ifirst + j;
                    for i in 0..len {
                        let k = (2 * j + i + pos) % nmidh + ifirst;
                        for l in 0..nc {
                            x[k][l] += tmp[jj][l] * g[i] + tmp[jj + nmid][l] * h[i];
                        }
                    }
                }
            }
        }
        for i in 0..n {
            for l in 0..nc {
                tmp[i][l] = x[i][l];
                x[i][l] = 0.0;
            }
        }
        let pos = position(n, filter.ishift);
        for j in 0..nhalf {
            for i in 0..len {
                let k = (2 * j + i + pos) % n;
                for l in 0..nc {
                    x[k][l] += tmp[j][l] * g[i] + tmp[j + nhalf][l] * h[i];
                }
            }
        }
    }
    1
}

fn trans_rows(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    filter: &WaveFilter,
    pow: usize,
    stage: usize,
    typ: i32,
) -> i32 {
    if stage > pow {
        return 0;
    }
    let n = 1 << pow;
    let nc = x[0].len();
    if stage == 0 {
        for i in 0..n {
            for l in 0..nc {
                if typ == 0 {
                    y[i][l] = x[i][l];
                } else {
                    x[i][l] = y[i][l];
                }
            }
        }
        return 1;
    }
    let len = filter.len;
    let h = &filter.filterh;
    let g = &filter.filterg;
    let nhalf = n >> 1;
    let mut tmp = vec![vec![0.0f32; nc]; n];
    if typ == 0 {
        let pos = position(n, filter.ishift);
        for j in 0..nhalf {
            for l in 0..nc {
                tmp[j][l] = 0.0;
                y[j][l] = 0.0;
            }
            for i in (1..len).step_by(2) {
                let k = (i + 2 * j - 1 + pos) % n;
                for l in 0..nc {
                    tmp[j][l] += x[k][l] * h[i - 1] + x[k + 1][l] * h[i];
                    y[j][l] += x[k][l] * g[i - 1] + x[k + 1][l] * g[i];
                }
            }
        }
        for ilevel in 1..stage {
            let nmid = 1 << (pow - ilevel);
            let nmidh = nmid >> 1;
            let ifirst = n - (nmid << 1);
            let ofirst = ifirst + nmid;
            let pos = position(nmid, filter.ishift);
            for j in 0..nmidh {
                let jj = ofirst + j;
                for l in 0..nc {
                    tmp[jj][l] = 0.0;
                    y[jj][l] = 0.0;
                }
                for i in (1..len).step_by(2) {
                    let k = (i + 2 * j - 1 + pos) % nmid + ifirst;
                    for l in 0..nc {
                        tmp[jj][l] += tmp[k][l] * h[i - 1] + tmp[k + 1][l] * h[i];
                        y[jj][l] += tmp[k][l] * g[i - 1] + tmp[k + 1][l] * g[i];
                    }
                }
            }
        }
        let nmid = 1 << (pow - stage);
        let ifirst = n - (nmid << 1);
        let ofirst = n - nmid;
        for i in 0..nmid {
            for l in 0..nc {
                y[ofirst + i][l] = tmp[ifirst + i][l];
            }
        }
    } else {
        for i in 0..n {
            for l in 0..nc {
                x[i][l] = y[i][l];
            }
        }
        for ilevel in (1..=stage).rev() {
            let nmid = 1 << (pow - ilevel);
            let nmidh = nmid << 1;
            let pos = position(nmidh, filter.ishift);
            let ifirst = n - nmidh;
            for i in 0..nmidh {
                for l in 0..nc {
                    tmp[ifirst + i][l] = x[ifirst + i][l];
                    x[ifirst + i][l] = 0.0;
                }
            }
            for j in 0..nmid {
                let jj = ifirst + j;
                for i in 0..len {
                    let k = (2 * j + i + pos) % nmidh + ifirst;
                    for l in 0..nc {
                        x[k][l] += tmp[jj][l] * g[i] + tmp[jj + nmid][l] * h[i];
                    }
                }
            }
        }
    }
    1
}

type RowKernel = fn(&mut [Vec<f32>], &mut [Vec<f32>], &WaveFilter, usize, usize, i32) -> i32;

fn one(
    x: &mut [f32],
    y: &mut [f32],
    filter: &WaveFilter,
    pow: usize,
    stage: usize,
    typ: i32,
    kernel: RowKernel,
) -> i32 {
    if stage > pow {
        return 0;
    }
    let n = 1 << pow;
    let mut xx = Vec::with_capacity(n);
    let mut yy = Vec::with_capacity(n);
    for i in 0..n {
        xx.push(vec![x[i]]);
        yy.push(vec![y[i]]);
    }
    let flag = kernel(&mut xx, &mut yy, filter, pow, stage, typ);
    for i in 0..n {
        x[i] = xx[i][0];
        y[i] = yy[i][0];
    }
    flag
}

/// C prototype: `int wavePack_1(float *x, float *y, waveFilter *filter, int twopow, int stage, int type);`
/// Transform length must be a power of two; `stage > twopow` returns zero.
pub fn wavePack_1(
    x: &mut [f32],
    y: &mut [f32],
    filter: &WaveFilter,
    twopow: usize,
    stage: usize,
    typ: i32,
) -> i32 {
    one(x, y, filter, twopow, stage, typ, pack_rows)
}

/// C prototype: `int waveTrans_1(float *x, float *y, waveFilter *filter, int twopow, int stage, int type);`
pub fn waveTrans_1(
    x: &mut [f32],
    y: &mut [f32],
    filter: &WaveFilter,
    twopow: usize,
    stage: usize,
    typ: i32,
) -> i32 {
    one(x, y, filter, twopow, stage, typ, trans_rows)
}

fn two(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    f1: &WaveFilter,
    f2: &WaveFilter,
    pow1: usize,
    pow2: usize,
    stage1: usize,
    stage2: usize,
    typ: i32,
    kernel: RowKernel,
) -> i32 {
    if stage1 > pow1 || stage2 > pow2 {
        return 0;
    }
    let n1 = 1 << pow1;
    let n2 = 1 << pow2;
    let mut z = vec![vec![0.0; n1]; n2];
    let mut flag = 1;
    if typ == 0 {
        for i in 0..n2 {
            if one(&mut x[i], &mut z[i], f1, pow1, stage1, typ, kernel) == 0 {
                flag = 0;
            }
        }
        if kernel(&mut z, y, f2, pow2, stage2, typ) == 0 {
            flag = 0;
        }
    } else {
        for i in 0..n2 {
            if one(&mut z[i], &mut y[i], f1, pow1, stage1, typ, kernel) == 0 {
                flag = 0;
            }
        }
        if kernel(x, &mut z, f2, pow2, stage2, typ) == 0 {
            flag = 0;
        }
    }
    flag
}

/// C prototype: `int wavePack_2(float **x, float **y, waveFilter *filter1, waveFilter *filter2, int twopow1, int twopow2, int stage1, int stage2, int type);`
pub fn wavePack_2(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    f1: &WaveFilter,
    f2: &WaveFilter,
    twopow1: usize,
    twopow2: usize,
    stage1: usize,
    stage2: usize,
    typ: i32,
) -> i32 {
    two(
        x, y, f1, f2, twopow1, twopow2, stage1, stage2, typ, pack_rows,
    )
}

/// C prototype: `int waveTrans_2(float **x, float **y, waveFilter *filter, int twopow1, int twopow2, int stage1, int stage2, int type);`
/// C returns EXIT_SUCCESS (zero), including for a successful transform.
pub fn waveTrans_2(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    filter: &WaveFilter,
    twopow1: usize,
    twopow2: usize,
    stage1: usize,
    stage2: usize,
    typ: i32,
) -> i32 {
    two(
        x, y, filter, filter, twopow1, twopow2, stage1, stage2, typ, trans_rows,
    );
    0
}
