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

use crate::util::{fieldvec_to_ringvec, hadamard, pad_matrtixrows_to, ConstraintSystemExt};

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

pub fn r1cs_to_ccs<F: PrimeField, R: Ring>(
    cs: ConstraintSystemRef<F>,
    z: &[R],
    l: usize,
    x_len: usize,
    wit_len: usize,
    W: usize,
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
    // let c_z = hadamard(a.clone(), b.clone(), c.clone(), z, new_r1cs_rows); // check

    // println!("rows padded!");
    let mut ccs = CCS {
        m: W,
        n: x_len + wit_len + 1,
        l: 1,
        t: 4,
        q: 2,
        d: 3,
        s: log2(W) as usize,
        s_prime: x_len + wit_len + 1,
        // M: vec![a, b, c, d_mat],
        // S: vec![vec![0, 1, 2], vec![3]],
        M: vec![a, b, c],
        S: vec![vec![0, 1], vec![2]],
        c: vec![R::one(), R::one().neg()],
    };
    let len = std::cmp::max((ccs.n - ccs.l - 1) * l, ccs.m).next_power_of_two();
    ccs.pad_rows_to(len);
    ccs
}

fn ret_ccs<
    const C: usize, // rows
    // const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
    F: PrimeField,
>(
    cs: ConstraintSystemRef<F>,
    x_r1cs: Vec<F>,
    w_r1cs: Vec<F>,
    wit_len: usize,
    w: usize
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, R>,
) {
    let mut rng = ark_std::test_rng();
    let x_ccs: Vec<R> = fieldvec_to_ringvec(&x_r1cs);
    let w_ccs: Vec<R> = fieldvec_to_ringvec(&w_r1cs);
    let one = R::one();

    let mut z = vec![one];
        z.extend(&x_ccs);
        z.extend(&w_ccs);
    let ccs: CCS<R> = r1cs_to_ccs::<F, R>(cs.clone(), &z, P::L, x_r1cs.len(), wit_len, w);
    ccs.check_relation(&z).expect("R1CS invalid!");
    println!("CCS check passed!\n");
    let scheme: AjtaiCommitmentScheme<C, R> = AjtaiCommitmentScheme::rand(&mut rng, w);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, P>(&scheme, w).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

pub fn setup_environment<
    const C: usize,
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    // W: usize,
    CS: LatticefoldChallengeSet<RqNTT>,
    F: PrimeField,
>(
    cs: ConstraintSystemRef<F>,
    w: usize
) -> (
    LCCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, RqNTT>,
) {
    let (x_len, wit_len, x_r1cs, w_r1cs) = z_split(cs.clone());
    let (cm_i, wit, ccs, scheme) =
        ret_ccs::<C, DP, RqNTT, F>(cs.clone(), x_r1cs, w_r1cs.clone(), wit_len, w);
    
    // aggregate the witnesses into w_ccs
    let ring_w_ccs:Vec<RqNTT>  =  wit.w_ccs
                                    .iter()
                                    .map(|elem| RqNTT::from(elem.clone()))
                                    .collect();
    let wit_acc = Witness::from_w_ccs::<DP>(ring_w_ccs);
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
