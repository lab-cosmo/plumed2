from typing import Dict, List, Optional

import torch
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatensor.torch.atomistic import (
    MetatensorAtomisticModel,
    ModelCapabilities,
    ModelMetadata,
    ModelOutput,
    System,
)
from rascaline.torch import SoapPowerSpectrum


class SOAP_CV(torch.nn.Module):
    def __init__(self, species):
        super().__init__()

        self.neighbor_type_pairs = Labels(
            names=["neighbor_1_type", "neighbor_2_type"],
            values=torch.tensor(
                [[t1, t2] for t1 in species for t2 in species if t1 <= t2]
            ),
        )
        self.calculator = SoapPowerSpectrum(
            cutoff=0.4,
            max_angular=2,
            max_radial=3,
            radial_basis={"Gto": {}},
            cutoff_function={"ShiftedCosine": {"width": 0.05}},
            center_atom_weight=1.0,
            atomic_gaussian_width=0.4,
        )

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels],
    ) -> Dict[str, TensorMap]:

        if "plumed::cv" not in outputs:
            return {}

        output = outputs["plumed::cv"]

        device = torch.device("cpu")
        if len(systems) > 0:
            device = systems[0].positions.device

        soap = self.calculator(systems, selected_samples=selected_atoms)
        soap = soap.keys_to_samples("center_type")
        soap = soap.keys_to_properties(self.neighbor_type_pairs)

        if not output.per_atom:
            raise ValueError("per_atom=False is not supported")

        soap_block = soap.block()

        samples = soap_block.samples.remove("center_type")

        block = TensorBlock(
            values=soap_block.values,
            samples=samples,
            components=[],
            properties=soap_block.properties,
        )
        cv = TensorMap(
            keys=Labels("_", torch.tensor([[0]], device=device)),
            blocks=[block],
        )

        return {"plumed::cv": cv}


cv = SOAP_CV(species=[8])
cv.eval()


capabilities = ModelCapabilities(
    outputs={
        "plumed::cv": ModelOutput(
            quantity="",
            unit="",
            per_atom=True,
            explicit_gradients=["postions"],
        )
    },
    interaction_range=0.4,
    supported_devices=["cpu", "mps", "cuda"],
    length_unit="nm",
    atomic_types=[8],
    dtype="float64",
)

metadata = ModelMetadata(
    name="Soap component retrieval",
    description="""
Retrieving all soap components for testing purposes
""",
    authors=["Gareth Tribello"],
    references={
        "implementation": ["ref to SOAP code"],
        "architecture": ["ref to SOAP"],
        "model": ["ref to paper"],
    },
)


model = MetatensorAtomisticModel(cv, metadata, capabilities)
model.export("soap_cv.pt", collect_extensions="extensions")
