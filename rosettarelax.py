import sys
import os
import pyrosetta

def Rosetta_relax(pdb_file):
    pyrosetta.init()
    from pyrosetta.rosetta.core.select import residue_selector as selections
    from pyrosetta import pose_from_pdb, create_score_function
    from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
    from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action
    from pyrosetta.rosetta.protocols.relax import FastRelax

    print(f'Rosetta processing {pdb_file} for Relax')

    pose = pose_from_pdb(pdb_file)
    scorefxn = create_score_function('ref2015')

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    tf.push_back(operation.PreventRepacking())

    all_residue_selector = selections.ResidueIndexSelector()
    if pose.num_chains() == 2:
        all_residue_selector.set_index_range(1, len(pose.chain_sequence(1)))
    else:
        all_residue_selector.set_index_range(1, len(pose.chain_sequence(1)) + len(pose.chain_sequence(2)))

    nbr_selector = selections.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(all_residue_selector)
    nbr_selector.set_include_focus_in_subset(True)
    subset_selector = nbr_selector

    prevent_repacking_rlt = operation.PreventRepackingRLT()
    prevent_subset_repacking = operation.OperateOnResidueSubset(
        prevent_repacking_rlt,
        subset_selector,
        flip_subset=True,
    )
    tf.push_back(prevent_subset_repacking)
    packer_task = tf.create_task_and_apply_taskoperations(pose)

    movemap = MoveMapFactory()
    movemap.add_bb_action(move_map_action.mm_enable, all_residue_selector)
    movemap.add_chi_action(move_map_action.mm_enable, subset_selector)
    mm = movemap.create_movemap_from_pose(pose)

    fastrelax = FastRelax()
    fastrelax.set_scorefxn(scorefxn)
    fastrelax.set_movemap(mm)
    fastrelax.set_task_factory(tf)
    fastrelax.apply(pose)

    pose.dump_pdb(pdb_file)
    print(f'Relaxed PDB saved to {pdb_file}')


def main():
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(__file__)} <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    if not os.path.isfile(pdb_file):
        print(f"Error: File '{pdb_file}' does not exist.")
        sys.exit(1)

    Rosetta_relax(pdb_file)

if __name__ == "__main__":
    main()
