use std::{fmt::Debug, marker::PhantomData, result, usize};

use ark_ff::PrimeField;
use ark_r1cs_std::{boolean::Boolean, eq::EqGadget};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError,
};
use ark_std::{log2, vec::Vec, UniformRand};
use arkworks_keccak::{
    constraints::{keccak_f_1600, pad101, KeccakMode},
    util::{bytes_to_bitvec, libary_step_sponge, vec_to_public_input},
};
use cyclotomic_rings::{challenge_set::LatticefoldChallengeSet, rings::SuitableRing};
use latticefold::{
    arith::{Arith, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::linearization::{LFLinearizationProver, LinearizationProver},
    transcript::poseidon::PoseidonTranscript,
};
use stark_rings_linalg::SparseMatrix;

use crate::util::{fieldvec_to_ringvec, ConstraintSystemExt};

pub struct KeccakCircuit<F: PrimeField, R: SuitableRing> {
    pub preimage: Vec<Boolean<F>>, // 512 bools
    pub expected: Vec<u8>,         // 32 bytes == 256 bits
    pub mode: KeccakMode,
    pub outputsize: usize, // binary output size
    pub r: usize,
    _phantom: PhantomData<R>,
}

pub struct StepCircuit<F: PrimeField, R: SuitableRing> {
    pub z_i: Vec<Boolean<F>>,
    pub expected_z_iplus1: Vec<Boolean<F>>,
    pub message_block: Option<Vec<Boolean<F>>>,
    pub flag: bool,
    pub r: usize,
    _phantom: PhantomData<R>,
}

impl<F: PrimeField, R: SuitableRing> ConstraintSynthesizer<F> for StepCircuit<F, R> {
    /// Create R1CS for z_iplus1 = F(z_i, w_i)
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        // let preimage: Vec<Boolean<F>> = vec_to_public_input(cs.clone(), "preimage", self.preimage)?;

        // let expected: Vec<Boolean<F>> = vec_to_public_input(cs.clone(), "expected", expected)?;
        // let result: Vec<Boolean<F>> =
        //     keccak_gadget(cs.clone(), &preimage, self.mode, self.outputsize)?;
        // Ensure constraints are generated without relying on n_steps
        // println!("Number of public inputs: {} + {}\n", preimage.len(), expected.len());
        // let (proof, result) =
        //     recursive_proof_result::<F, R>(cs.clone(), n_steps, preimage, self.mode, self.outputsize)?;
        let z_i: Vec<Boolean<F>> = vec_to_public_input(cs.clone(), "z_i", self.z_i)?;
        let expected: Vec<Boolean<F>> = self.expected_z_iplus1;
        let z_iplus1: Vec<Boolean<F>> = keccak_f_1600(cs.clone(), &z_i)?;

        for (o, e) in z_iplus1.iter().zip(expected.iter()) {
            o.enforce_equal(e)?;
        }

        Ok(())
    }
}

impl<F: PrimeField, R: SuitableRing> KeccakCircuit<F, R> {
    pub fn new(
        preimage: Vec<Boolean<F>>,
        expected: Vec<u8>,
        mode: KeccakMode,
        outputsize: usize,
        r: usize,
    ) -> Self {
        Self {
            preimage,
            expected,
            mode,
            outputsize,
            r,
            _phantom: PhantomData,
        }
    }
}

impl<F: PrimeField, R: SuitableRing> StepCircuit<F, R> {
    pub fn new(
        z_i: Vec<Boolean<F>>,
        expected_z_iplus1: Vec<Boolean<F>>,
        message_block: Option<Vec<Boolean<F>>>,
        flag: bool,
        r: usize,
    ) -> Self {
        Self {
            z_i,
            expected_z_iplus1,
            message_block,
            flag,
            r,
            _phantom: PhantomData,
        }
    }
}

pub trait InitCircuit {
    fn init_constraint_system<F: PrimeField>() -> ConstraintSystemRef<F>;
}

pub struct CircuitInit;

impl InitCircuit for CircuitInit {
    fn init_constraint_system<F: PrimeField>() -> ConstraintSystemRef<F> {
        use ark_relations::r1cs::{ConstraintLayer, ConstraintSystem, TracingMode};
        use tracing_subscriber::{layer::SubscriberExt, Registry};

        let mut layer = ConstraintLayer::default();
        layer.mode = TracingMode::OnlyConstraints;
        let subscriber = Registry::default().with(layer);
        let _guard = tracing::subscriber::set_default(subscriber);

        let cs = ConstraintSystem::new_ref();
        cs
    }
}

pub fn z_split_lengths<F: PrimeField>(cs: ConstraintSystemRef<F>) -> (usize, usize) {
    let public_inputs = cs.ret_instance();
    let witnesses = cs.ret_witness();
    assert_eq!(public_inputs[0], F::one());
    println!(
        "Witnesses: {} + {} = {}\n",
        public_inputs.len(),
        witnesses.len(),
        public_inputs.len() + witnesses.len()
    );

    (public_inputs.len() - 1, witnesses.len())
}

pub fn r1cs_to_ccs<F: PrimeField, R: SuitableRing>(
    cs: ConstraintSystemRef<F>,
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

pub fn ret_ccs<
    const C: usize, // rows
    // const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
    F: PrimeField,
>(
    cs: ConstraintSystemRef<F>,
    w: usize,
    scheme: AjtaiCommitmentScheme<C, R>
) -> (CCCS<C, R>, Witness<R>, CCS<R>) {
    let x_ccs: Vec<R> = fieldvec_to_ringvec(&cs.ret_instance()[1..]);
    let w_ccs: Vec<R> = fieldvec_to_ringvec(&cs.ret_witness());

    // z = 1 || x || w
    // println!("Length of F::one in ring vec: {}",fieldvec_to_ringvec::<R, F>(&[F::from(1u128)]).len());
    let mut z: Vec<R> = Vec::with_capacity(x_ccs.len() + w_ccs.len() + 1);
    z.extend_from_slice(&vec![fieldvec_to_ringvec(&[F::from(1u128)])[0]]);
    z.extend_from_slice(&x_ccs);
    z.extend_from_slice(&w_ccs);

    let ccs: CCS<R> = r1cs_to_ccs::<F, R>(cs.clone(), P::L, x_ccs.len(), w_ccs.len(), w);
    ccs.check_relation(&z).expect("R1CS invalid!");
    println!("CCS<R> check passed!\n");

    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, P>(&scheme, w).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs)
}

pub fn ret_linearized_inst_wit_pair<
    F: PrimeField,
    R: SuitableRing,
    const C: usize,
    DP: DecompositionParams,
    CS: LatticefoldChallengeSet<R>,
>(
    wit: Witness<R>,
    cm_i: CCCS<C, R>,
    ccs: CCS<R>,
) -> Result<(LCCCS<C, R>, Witness<R>), SynthesisError> {
    let agrt_w_ccs: Vec<R> = wit.w_ccs.clone();
    let wit_acc = Witness::from_w_ccs::<DP>(agrt_w_ccs);
    let mut transcript = PoseidonTranscript::<R, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<R, CS>>::prove::<C, F>(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .expect("Failed to generate linearization proof");
    Ok((acc, wit_acc))
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
    w: usize,
) -> (
    LCCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, RqNTT>,
) {
    let mut rng = ark_std::rand::thread_rng();
    let scheme = AjtaiCommitmentScheme::rand(&mut rng, w);
    let (cm_i, wit, ccs) = ret_ccs::<C, DP, RqNTT, F>(cs.clone(), w, scheme.clone());

    // reference to CCS witness for cm_i
    let agrt_w_ccs: Vec<RqNTT> = wit.w_ccs.clone();
    let wit_acc = Witness::from_w_ccs::<DP>(agrt_w_ccs);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove::<C, F>(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .expect("Failed to generate linearization proof");

    (acc, wit_acc, cm_i, wit, ccs, scheme)
}
