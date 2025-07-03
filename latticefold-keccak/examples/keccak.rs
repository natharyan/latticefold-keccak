use std::{fmt::Debug, result, time::Instant};

use ark_bls12_381::Fr;
use ark_ff::PrimeField;
use ark_r1cs_std::prelude::Boolean;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{rand::Rng, vec::Vec};
use arkworks_keccak::{
    constraints::{keccak_f_1600, pad101, ret_r, split_to_blocks, truncate, KeccakMode},
    util::{bytes_to_bitvec, libary_step_sponge, shake_128},
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing},
};
use latticefold::{
    arith::{error::CSError, Arith, Witness, CCCS, CCS, LCCCS}, commitment::AjtaiCommitmentScheme, nifs::{LFProof, NIFSProver, NIFSVerifier}, transcript::{poseidon::PoseidonTranscript, TranscriptWithShortChallenges}
};
use latticefold_keccak::{
    circuit::{
        ret_ccs, ret_linearized_inst_wit_pair, setup_environment, z_split_lengths, CircuitInit,
        InitCircuit, KeccakCircuit, StepCircuit,
    },
    util::{ConstraintSystemExt},
};

include!("../src/constants.rs");

type RqNTT = GoldilocksRingNTT;
type CS = GoldilocksChallengeSet;
type T = PoseidonTranscript<RqNTT, CS>;

fn num_steps(num_blocks: usize, d: usize, r: usize) -> usize {
    num_blocks + (d / r + 1)
}

struct FoldingWitness<F: PrimeField> {
    message: Vec<Boolean<F>>,
    flag: bool,
}
pub struct RecursiveSnark<F, NTT, CS, DP, PsT>
where
    F: PrimeField,
    NTT: SuitableRing,
    CS: LatticefoldChallengeSet<NTT>,
    DP: DecompositionParams,
    PsT: TranscriptWithShortChallenges<NTT>
{
    _marker: std::marker::PhantomData<(F, NTT, CS, DP, PsT)>,
}

impl<
        F: PrimeField,
        NTT: SuitableRing,
        CS: LatticefoldChallengeSet<NTT>,
        DP: DecompositionParams,
        PsT: TranscriptWithShortChallenges<NTT>
    > RecursiveSnark<F, NTT, CS, DP, PsT>
{
    /// Returns a CS, cm_i, wit_i for step i and updates z_i and message_block
    fn prove_step(
        mut step_circuit: StepCircuit<F, NTT>,
        message_i: Option<Vec<Boolean<F>>>,
        flag_next: bool,
        mut folded_acc_i: Option<LCCCS<C, NTT>>,
        mut folded_witness_i: Option<Witness<NTT>>,
        mut output: Vec<Boolean<F>>,
        i: usize
    ) -> Result<(StepCircuit<F, NTT>, LCCCS<C, NTT>, Witness<NTT>, LFProof<C, NTT>, Vec<Boolean<F>>), SynthesisError> {

        let cs: ConstraintSystemRef<F> = CircuitInit::init_constraint_system::<F>();

        let flag_f = if step_circuit.flag {F::one()} else {F::zero()};
        let flag_next_f = if flag_next {F::one()} else {F::zero()};
        // Ensure that if the current flag == 1, then flag_next == 1. If flag_current == 0, then flag_next can be 0 or 1
        assert!(flag_f * (F::one() - flag_next_f) == F::zero());
        
        step_circuit.message_block = message_i;
        // Compute expected_z_iplus1 for the current step
        let expected_z_iplus1 = libary_step_sponge(
            step_circuit.z_i.clone(),
            step_circuit.message_block.clone(),
            step_circuit.r,
            Boolean::Constant(step_circuit.flag),
        )?;
        step_circuit.expected_z_iplus1 = expected_z_iplus1;

        // Absorption phase, else Squeezing phase
        if step_circuit.flag == false {
            let mut z_i = step_circuit.z_i.clone();
            if let Some(message_block) = &step_circuit.message_block {
                for i in 0..step_circuit.r {
                    z_i[i] = Boolean::xor(&z_i[i], &message_block[i])?;
                }
            } else {
                return Err(SynthesisError::AssignmentMissing);
            }
            step_circuit.z_i = z_i;
        } 

        let sc = StepCircuit::<F, NTT>::new(
            step_circuit.z_i.clone(),
            step_circuit.expected_z_iplus1.clone(),
            step_circuit.message_block.clone(),
            step_circuit.flag,
            step_circuit.r,
        );

        sc.generate_constraints(cs.clone())?;

        // check if the constraint system is satisfied
        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied{
            println!("Unsatisfied constraint: {:?}\n", cs.which_is_unsatisfied());
        }
        assert!(is_satisfied);

        let mut prover_transcript = PoseidonTranscript::<NTT, CS>::default();
        let w: usize;

        let wit_len = cs.ret_witness().len();
        w = wit_len * DP::L;
        let mut rng = ark_std::rand::thread_rng();
        let scheme_i = AjtaiCommitmentScheme::rand(&mut rng, w);

        // if i == 0, then initialize folded_acc and folded_witness
        if i == 0{
            let (cm_init, wit_init, ccs) = ret_ccs::<C, DP, NTT, F>(cs.clone(), w, scheme_i.clone());
            let (acc, witness) = ret_linearized_inst_wit_pair::<F, NTT, C, DP, CS>(wit_init, cm_init, ccs)?;
            folded_acc_i = Some(acc);
            folded_witness_i = Some(witness);
        }
        
        // Generate witness and commitment for the current step
        let (cm_i, wit_i, ccs) = ret_ccs::<C, DP, NTT, F>(cs.clone(), w, scheme_i.clone());

        let (folded_acc_iplus1, folded_witness_iplus1, proof) = NIFSProver::<C, NTT, DP, PsT>::prove::<Fr>(
            &folded_acc_i.as_ref().expect("folded_acc unitialized"),
            &folded_witness_i.as_ref().expect("folded_witness unitialized"),
            &cm_i,
            &wit_i,
            &mut prover_transcript,
            &ccs,
            &scheme_i,
            w,
        )
        .unwrap();

        
        if step_circuit.flag == true {
            output.extend(truncate(&step_circuit.expected_z_iplus1, step_circuit.r)?);
        }

        // Update the circuit
        step_circuit.z_i = step_circuit.expected_z_iplus1.clone();

        // Ensure output is correct
        if !output.is_empty() {
            let last_r_bits = &output[output.len() - step_circuit.r..];
            assert_eq!(last_r_bits, &truncate(&step_circuit.expected_z_iplus1, step_circuit.r)?);
        }

        Ok((step_circuit, folded_acc_iplus1, folded_witness_iplus1, proof, output))
    }

    // call recursive_proof_result; for each step create a separate

    /// public inputs: z_0, n_steps, r
    /// witnesses - message blocks, flag, (folded_acc, folded_witness)
    fn recursive_proof_result(
        circuit: KeccakCircuit<F, NTT>,
        n_steps: usize,
        z_0: Vec<Boolean<F>>,
        output: Vec<Boolean<F>>,
        message: Vec<Boolean<F>>,
        mode: KeccakMode,
        outputsize: usize,
    ) -> Result<(LFProof<C, NTT>, Vec<Boolean<F>>), SynthesisError> {
        // TODO: init flag as false, don't allow it to be changed to true once changed to false.
        // assert and of not flag and message_block is None

        let m_blocks = split_to_blocks(&message, circuit.r)?;
        let mut step_circuit_i =
            StepCircuit::<F, NTT>::new(z_0, vec![], None, false, circuit.r);

        let mut flag_next = false;

        let mut proof: Option<LFProof<C, NTT>> = None;

        let mut folded_acc: Option<LCCCS<C, NTT>> = None;
        let mut folded_witness: Option<Witness<NTT>> = None;
        
        for i in 0..n_steps {
            if (i + 1) % m_blocks.len() == 0 {
                flag_next = true;
            }
            let start_step = Instant::now();
            // Call step_function for the current step
            let (step_circuit_iplus1, folded_acc_iplus1, folded_witness_iplus1, folding_proof, output) = Self::prove_step(step_circuit_i, Some(m_blocks[i].clone()), flag_next, folded_acc, folded_witness, output.clone(), i)?;
            proof = Some(folding_proof);
            folded_acc = Some(folded_acc_iplus1);
            folded_witness = Some(folded_witness_iplus1);
            step_circuit_i = step_circuit_iplus1;
        }

        Ok((proof.expect("proof must exist"), output))
    }
}

fn main() {
    let mut rng = ark_std::rand::thread_rng();
    let preimage_length_bytes = rng.gen_range(1..=10);
    let preimage: Vec<u8> = (0..preimage_length_bytes).map(|_| rng.r#gen()).collect();
    let d: usize = rng.gen_range(50..=100);
    println!("input length: {} bytes", preimage.len());
    println!("d: {} bits", d);
    let expected = shake_128(&preimage, d / 8); // change
    let preimage = bytes_to_bitvec::<Fr>(&preimage);

    
    let mode: KeccakMode = KeccakMode::Shake128;
    let r = ret_r(mode);
    let padded_preimage = pad101(&preimage, mode).unwrap();
    let num_blocks = split_to_blocks(&padded_preimage, r).iter().len();

    let n_steps = num_steps(num_blocks, d, r);

    // TODO: change d as per KeccakMode
    let circuit = KeccakCircuit::<Fr, RqNTT>::new(
        preimage.clone(),
        expected.to_vec(),
        KeccakMode::Shake128,
        d,
        r,
    );

    let z_0: Vec<Boolean<Fr>> = vec![Boolean::Constant(false); 1600]; // b = 1600
    let mut z_n: Vec<Boolean<Fr>> = vec![];

    let recursive_result = RecursiveSnark::<Fr, RqNTT, CS, GoldilocksExampleDP, T>::recursive_proof_result(circuit, n_steps, preimage, z_0, z_n, mode, d).unwrap();
    let proof = recursive_result.0;
    let z_n = recursive_result.1;

    
    // ///////////// TODO: Change from here: /////////////////////////////
    // circuit.generate_constraints();
    // let is_satisfied = cs.is_satisfied().unwrap();
    // assert!(is_satisfied);

    // let (x_len, w_len) = z_split_lengths(cs.clone());
    // assert_eq!(x_len, preimage_length_bytes * 8 + expected.len() * 8);

    // println!("Setting up example environment...");

    // println!("Decomposition parameters:");

    // println!("\tB: {}", GoldilocksExampleDP::B);
    // println!("\tL: {}", GoldilocksExampleDP::L);
    // println!("\tB_SMALL: {}", GoldilocksExampleDP::B_SMALL);
    // println!("\tK: {}", GoldilocksExampleDP::K);

    // let (acc, wit_acc, cm_i, wit_i, ccs, scheme) =
    //     setup_environment::<C, RqNTT, GoldilocksExampleDP, CS, Fr>(
    //         cs,
    //         w_len * GoldilocksExampleDP::L,
    //     );

    // let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    // let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    // println!("Generating proof...");
    // let start = Instant::now();
    // // TODO: convert to recursive proof by using folded_acc and folded_wit at each step instead of acc and wit_acc
    // let (folded_acc, folded_wit, proof) =
    //     NIFSProver::<C, RqNTT, GoldilocksExampleDP, T>::prove::<Fr>(
    //         &acc,
    //         &wit_acc,
    //         &cm_i,
    //         &wit_i,
    //         &mut prover_transcript,
    //         &ccs,
    //         &scheme,
    //         w_len * GoldilocksExampleDP::L,
    //     )
    //     .unwrap();
    // let duration = start.elapsed();
    // println!("Proof generated in {:?}", duration);

    // let mut serialized_proof = Vec::new();

    // println!("Serializing proof (with compression)...");
    // proof
    //     .serialize_with_mode(&mut serialized_proof, Compress::Yes)
    //     .unwrap();
    // let compressed_size = serialized_proof.len();
    // println!(
    //     "Proof size (with compression) size: {}",
    //     humansize::format_size(compressed_size, humansize::BINARY)
    // );

    // println!("Serializing proof (without compression)...");
    // proof
    //     .serialize_with_mode(&mut serialized_proof, Compress::No)
    //     .unwrap();
    // let uncompressed_size = serialized_proof.len();
    // println!(
    //     "Proof (without compression) size: {}\n",
    //     humansize::format_size(uncompressed_size, humansize::BINARY)
    // );

    // println!("Verifying proof");
    // let start = Instant::now();
    // // assert_eq!(folded_acc, acc, "folded_acc != acc");
    // NIFSVerifier::<C, RqNTT, GoldilocksExampleDP, T>::verify::<Fr>(
    //     &acc,
    //     &cm_i,
    //     &proof,
    //     &mut verifier_transcript,
    //     &ccs,
    // )
    // .unwrap();
    // let duration = start.elapsed();
    // println!("Proof verified in {:?}", duration);
}
