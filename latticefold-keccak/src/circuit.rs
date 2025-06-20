use std::{fmt::Debug, time::Instant, usize};

use ark_bls12_381::Fr;
use ark_ff::PrimeField;
use ark_relations::r1cs::{self, ConstraintSystemRef, Field};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{log2, rand::Rng, vec::Vec, UniformRand};
use arkworks_keccak::{
    constraints::{KeccakCircuit, KeccakMode},
    util::{bytes_to_bitvec, keccak256, sha3_256, shake_128, shake_256},
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{BabyBearChallengeSet, BabyBearRingNTT, SuitableRing},
};
use latticefold::{
    arith::{Arith, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};
use stark_rings::Ring;
use stark_rings_linalg::SparseMatrix;

use crate::util::{fieldvec_to_ringvec, hadamard_ret_d, pad_matrtixrows_to, ConstraintSystemExt};

pub fn z_split<F: PrimeField>(cs: ConstraintSystemRef<F>) -> (usize, usize, Vec<F>, Vec<F>) {
    let public_inputs = cs.ret_instance();
    let witnesses = cs.ret_witness();
    assert_eq!(public_inputs[0], F::one());
    println!("Witnesses: {}\n", public_inputs.len() + witnesses.len());

    (
        public_inputs.len() - 1,
        witnesses.len(),
        public_inputs[1..].to_vec(),
        witnesses,
    )
}

pub fn r1cs_to_ccs<F: PrimeField, R: Ring, const W: usize>(
    cs: ConstraintSystemRef<F>,
    z: &[R],
    l: usize,
    x_len: usize,
    wit_len: usize,
) -> CCS<R> {
    // D is the haddamard product of the matrices
    let (mut a, mut b, mut c): (SparseMatrix<R>, SparseMatrix<R>, SparseMatrix<R>) =
        cs.get_r1cs_matrices();
    let r1cs_rows = a.n_rows;
    let new_r1cs_rows = if l == 1 && (wit_len > 0 && (wit_len & (wit_len - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows
    };
    let d_mat = hadamard_ret_d(a.clone(), b.clone(), c.clone(), z, new_r1cs_rows);

    a = pad_matrtixrows_to(a.clone(), c.n_cols);
    b = pad_matrtixrows_to(b.clone(), c.n_cols);
    let mut ccs = CCS {
        m: W,
        n: z.len(),
        l: 1,
        t: 4,
        q: 2,
        d: 3,
        s: log2(W) as usize,
        s_prime: z.len(),
        M: vec![a, b, c, d_mat],
        S: vec![vec![0, 1, 2], vec![3]],
        c: vec![R::one(), R::one().neg()],
    };
    let len = std::cmp::max((ccs.n - ccs.l - 1) * l, ccs.m).next_power_of_two();
    ccs.pad_rows_to(len);
    ccs
}

fn ret_ccs<
    const C: usize, // rows
    const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
    F: PrimeField,
>(
    cs: ConstraintSystemRef<F>,
    x_r1cs: Vec<F>,
    w_r1cs: Vec<F>,
    wit_len: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, W, R>,
) {
    let mut rng = ark_std::test_rng();
    let x_ccs: Vec<R> = fieldvec_to_ringvec(&x_r1cs);
    let w_ccs: Vec<R> = fieldvec_to_ringvec(&w_r1cs);
    let one = R::one();

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    let ccs: CCS<R> = r1cs_to_ccs::<F, R, W>(cs.clone(), &z, P::L, x_r1cs.len(), wit_len);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

pub fn setup_environment<
    const C: usize,
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const W: usize,
    CS: LatticefoldChallengeSet<RqNTT>,
    F: PrimeField,
>(
    cs: ConstraintSystemRef<F>,
) -> (
    LCCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let (x_len, wit_len, x_r1cs, w_r1cs) = z_split(cs.clone());
    let (cm_i, wit, ccs, scheme) =
        ret_ccs::<C, W, DP, RqNTT, F>(cs.clone(), x_r1cs, w_r1cs.clone(), wit_len);
    let ring_w_css: Vec<RqNTT> = w_r1cs
        .clone() // TODO: match with the linearization protocol
        .into_iter()
        .map(|bit| {
            // expect boolean witnesses from r1cs-std
            if bit == F::zero() {
                RqNTT::from(0 as u64)
            } else {
                RqNTT::from(1 as u64)
            }
        })
        .collect();
    let wit_acc = Witness::from_w_ccs::<DP>(ring_w_css);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .expect("Failed to generate linearization proof");

    (acc, wit_acc, cm_i, wit, ccs, scheme)
}
