#ifdef __PLUMED_HAS_RASCALINE
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

#include <rascaline.hpp>


namespace PLMD { namespace rascaline {

/// Implementation of the `rascaline::System` interface containing data from
/// Plumed
class PlumedSystem final: public ::rascaline::System {
public:
    uintptr_t size() const override {
        return species_.size();
    }

    const int32_t* species() const override {
        return species_.data();
    }

    const double* positions() const override {
        return positions_.data();
    }

    CellMatrix cell() const override {
        return cell_;
    }

    void compute_neighbors(double cutoff) override {
        throw std::runtime_error(
            "can not use PlumedSystem to compute neighbors, set use_native_system to true"
        );
    }

    const std::vector<rascal_pair_t>& pairs() const override {
        throw std::runtime_error(
            "can not use PlumedSystem to compute neighbors, set use_native_system to true"
        );
    }

    const std::vector<rascal_pair_t>& pairs_containing(uintptr_t center) const override {
        throw std::runtime_error(
            "can not use PlumedSystem to compute neighbors, set use_native_system to true"
        );
    }

    std::vector<int32_t> species_;
    std::vector<double> positions_;
    CellMatrix cell_;
};

class RascalinePlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit RascalinePlumedAction(const ActionOptions&);

    // active methods:
    unsigned getNumberOfDerivatives() const override;
    void calculate() override;
    void apply() override;

private:
    void update_system();


    std::string calculator_name_;
    std::string hyper_parameters_;
    ::rascaline::Calculator calculator_;

    PlumedSystem system_;

    equistore::TensorMap descriptor_;
};

PLUMED_REGISTER_ACTION(RascalinePlumedAction, "RASCALINE")

void RascalinePlumedAction::registerKeywords(Keywords& keys) {
    Action::registerKeywords(keys);
    ActionAtomistic::registerKeywords(keys);
    ActionWithValue::registerKeywords(keys);

    keys.add("numbered", "SPECIES", "the atoms in each species type");
    keys.reset_style("SPECIES","atoms");

    keys.add("compulsory", "CALCULATOR", "name of the rascaline calculator to use");
    keys.add("compulsory", "PARAMETERS", "JSON string containing the hyper-parameters for the calculator");
}

std::string read_key(Action& action, const std::string& key) {
    std::string buffer;
    action.parse(key, buffer);
    return buffer;
}

RascalinePlumedAction::RascalinePlumedAction(const ActionOptions& ao):
    Action(ao),
    ActionAtomistic(ao),
    ActionWithValue(ao),
    calculator_name_(read_key(*this, "CALCULATOR")),
    hyper_parameters_(read_key(*this, "PARAMETERS")),
    calculator_(calculator_name_, hyper_parameters_),
    system_(),
    descriptor_(equistore::Labels({"_"}, {}), {}) // initialize with an empty TensorMap
{
    // Now read in atoms that we are using
    std::vector<AtomNumber> all_atoms;
    this->parseAtomList("SPECIES", all_atoms);

    if (all_atoms.empty()) {
        // parse multiple species lists
        std::vector<AtomNumber> current_atoms_list;
        int current_species = 0;

        while (true) {
            current_species += 1;
            current_atoms_list.clear();

            this->parseAtomList("SPECIES", current_species, current_atoms_list);
            if (current_atoms_list.empty()) {
                break;
            }

            log.printf("  Species %d includes atoms : ", current_species);
            for (const auto& atom: current_atoms_list) {
                log.printf("%d ", atom.serial());
                system_.species_.push_back(current_species);
                all_atoms.push_back(atom);
            }
            log.printf("\n");
        }
    } else {
        // we only have one species
        log.printf("  Species 1 includes atoms : ");
        for (const auto& atom: all_atoms) {
            log.printf("%d ", atom.serial());
            system_.species_.push_back(1);
        }
        log.printf("\n");
    }
    system_.positions_.resize(3 * all_atoms.size());

    // Request the atoms from plumed
    this->requestAtoms(all_atoms);

    // declare the value we will use to store the representation
    // NOTE: we declare it to be N x 1, and we will resize it later
    // to the actual size for the current representation
    this->addValue({static_cast<unsigned>(all_atoms.size()), 1});
    this->setNotPeriodic();
    this->getPntrToComponent(0)->alwaysStoreValues();

    this->checkRead();
}

unsigned RascalinePlumedAction::getNumberOfDerivatives() const {
    // forces (3 x N values) + stress (9 values)
    return 3 * this->getNumberOfAtoms() + 9;
}

void RascalinePlumedAction::update_system() {
    // Conversion factor from PLUMED length units to angstroms
    double plumed_to_angstroms = 1.0;
    if (!plumed.usingNaturalUnits()) {
        plumed_to_angstroms = 10 * plumed.getUnits().getLength();
    }

    for (size_t i=0; i<this->getNumberOfAtoms(); i++) {
        auto position = this->getPosition(i);
        system_.positions_[3 * i + 0] = plumed_to_angstroms * position[0];
        system_.positions_[3 * i + 1] = plumed_to_angstroms * position[1];
        system_.positions_[3 * i + 2] = plumed_to_angstroms * position[2];
    }

    // Set the box and the pbc from the cell
    auto box = this->getBox();
    for(size_t i=0; i<3; i++) {
        for(size_t j=0; j<3; j++) {
            system_.cell_[i][j]= plumed_to_angstroms * box(i, j);
        }
    }
}


void RascalinePlumedAction::calculate() {
    this->update_system();

    auto options = ::rascaline::CalculationOptions();
    options.use_native_system = true;

    // only compute derivatives if needed
    if (!this->doNotCalculateDerivatives()) {
        options.gradients = {"positions", "cell"};
    }

    auto systems = std::vector<::rascaline::System*>{&system_};
    descriptor_ = calculator_.compute(systems, options);

    if (calculator_name_ == "spherical_expansion") {
        descriptor_ = descriptor_.keys_to_properties("species_neighbor");
        descriptor_ = descriptor_.components_to_properties("spherical_harmonics_m");
        descriptor_ = descriptor_.keys_to_properties("spherical_harmonics_l");
        descriptor_ = descriptor_.keys_to_samples("species_center", /*sort_samples*/ true);
    } else if (calculator_name_ == "soap_power_spectrum") {
        auto neighbors = std::vector<std::string>{"species_neighbor_1", "species_neighbor_2"};
        descriptor_ = descriptor_.keys_to_properties(neighbors);
        descriptor_ = descriptor_.keys_to_samples("species_center", /*sort_samples*/ true);
    } else {
        // TODO: add the other relevant representations here

        throw std::runtime_error(
            "unknown calculator " + calculator_name_ + " edit this file to add "
            "new calculators"
        );
    }

    // there should be a single block now
    assert(descriptor_.keys().count() == 1);
    auto block = descriptor_.block_by_id(0);

    auto data = block.values();
    auto shape = data.shape();

    auto value = this->getPntrToComponent(0);
    assert(value->getShape()[0] == shape[0]);

    // reshape the output `Value` if needed
    if (value->getShape()[1] != shape[1]) {
        value->setShape({this->getNumberOfAtoms(), static_cast<unsigned>(shape[1])});
    }

    for (size_t i=0; i<shape[0]; i++) {
        for (size_t j=0; j<shape[1]; ++j) {
            value->set(i * shape[1] + j, data(i, j));
        }
    }
}


void RascalinePlumedAction::apply() {
    throw "not implemented";
//   // Do nothing if no forces were added
//   if( !getPntrToOutput(0)->forcesWereAdded() ) return ;
//   // Check that plumed can handle the forces. The CV doesn't work if the box is too small
//   for(unsigned i=0; i<3; ++i) {
//       if( getBox()[i][i]<cutoff ) error("cannot calculate rascal derivatives correctly for small cells");
//   }
//   // Clear the forces from last time
//   std::fill(forcesToApply.begin(),forcesToApply.end(),0);

//   // Recalculate SOAPS if forces were added
//   structureToJson();
//   // Now lets setup the stuff for librascal
//   auto manager = make_structure_manager_stack<StructureManagerCenters, AdaptorNeighbourList, AdaptorCenterContribution, AdaptorStrict>( structure, adaptors );
//   // And compute everything with librascal
//   representation->compute(manager);

//   // Now lets get the gradients
//   auto && soap_vector_gradients{*manager->template get_property<PropGrad_t>(representation->get_gradient_name())};
//   math::Matrix_t gradients = soap_vector_gradients.get_features_gradient();
//   // This gets information on the neighbour lists from librascal
//   Eigen::Matrix<int, Eigen::Dynamic, 5> ninfo = manager->get_gradients_info();
//   // Determine the number of neighbours for each of the environments
//   std::fill(neigh.begin(),neigh.end(),0); for(unsigned i=0;i<ninfo.rows();++i) neigh[ ninfo(i,1) ]++;

//   // Apply the forces on the atoms
//   Value* outval=getPntrToOutput(0); Tensor vir; vir.zero();
//   double lunit = 1.0; if( !plumed.usingNaturalUnits() ) lunit = 10*plumed.getUnits().getLength();
//   std::vector<unsigned> shape( outval->getShape() ); unsigned base=0;
//   // Loop over environments (atoms)
//   for(unsigned i=0; i<shape[0]; ++i) {
//       // Loop over neighbours in ith environment
//       for(unsigned k=0; k<neigh[i]; ++k) {
//           // Get distance for calculation of virial contribution
//           Vector force, dist = pbcDistance( getPosition(i), getPosition(ninfo(base+k,2)) );
//           // Loop over features
//           for(unsigned j=0; j<shape[1];++j) {
//               // Get the force on the jth feature in the ith environment.
//               double ff = outval->getForce( i*shape[1] + j );
//               // Loop over x, y, z
//               for(unsigned n=0;n<3;++n) {
//                   force[n] = ff*gradients(3*(base+k)+n,j)*lunit;
//                   // This is the force to add to the central atom
//                   forcesToApply[ 3*i+n ] -= force[n];
//                   // This is the force to add to the neighbour atom.
//                   forcesToApply[ 3*ninfo(base+k,2) + n ] += force[n];
//               }
//               vir -= Tensor(force,dist);
//           }
//       }
//       // And this is some book keeping to make sure we are navigating gradients and info correctly
//       base = base + neigh[i];
//   }
//   unsigned vbase = 3*getNumberOfAtoms();
//   for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) forcesToApply[vbase + 3*i + j] = vir(i,j);
//   unsigned start=0; setForcesOnAtoms( forcesToApply, start );
}


} } // namespace PLMD; namespace rascaline


#endif
