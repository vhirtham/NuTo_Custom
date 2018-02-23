#include <boost/test/unit_test.hpp>

// NuTo includes
#include "mechanics/cell/SimpleAssembler.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/mesh/UnitMeshFem.h"

// NuToCustom includes
#include "integrands/MoistureTransport.h"

// Other includes
#include <boost/ptr_container/ptr_vector.hpp>

using namespace NuTo;
using namespace NuTo::Integrands;
using namespace std::placeholders;


// Custom moisture transport coefficients %%%%%%%

class DCw
{
    constexpr static double diffusionCoefficient = 1;

public:
    constexpr static double value(double WV, double RH, double WV_dt, double RH_dt)
    {
        return diffusionCoefficient;
    }
    constexpr static double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
    constexpr static double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
};


class DCg
{
    constexpr static double diffusionCoefficient = 1;

public:
    constexpr static double value(double WV, double RH, double WV_dt, double RH_dt)
    {
        return diffusionCoefficient;
    }
    constexpr static double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
    constexpr static double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
};


class MeC
{
    constexpr static double massExchangeCoefficient = 1;

public:
    constexpr static double value(double WV, double RH, double WV_dt, double RH_dt)
    {
        return massExchangeCoefficient;
    }
    constexpr static double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
    constexpr static double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
};

class WVEq
{
    constexpr static double coeff0 = 0.;
    constexpr static double coeff1 = 0.18;

public:
    constexpr static double value(double WV, double RH, double WV_dt, double RH_dt)
    {
        return coeff1 * RH + coeff0;
    }
    constexpr static double d_WaterVolumeFraction(double WV, double RH, double WV_dt, double RH_dt)
    {
        return 0.;
    }
    constexpr static double d_RelativeHumidity(double WV, double RH, double WV_dt, double RH_dt)
    {
        return coeff1;
    }
};

BOOST_AUTO_TEST_CASE(IntegrationTest)
{
    // Create mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constexpr int numElements = 2;
    MeshFem mesh = UnitMeshFem::CreateLines(numElements);

    DofType dofRH("relativeHumidity", 1);
    DofType dofWV("waterVolumeFraction", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLobatto(2));
    AddDofInterpolation(&mesh, dofRH, interpolation);
    AddDofInterpolation(&mesh, dofWV, interpolation);


    // Create constraints %%%%%%%%%%%%%%%%%%%%%%%
    Constraint::Constraints constraints;

    // Set nodal values %%%%%%%%%%%%%%%%%%%%%%%%%
    Group<NodeSimple> groupNodesRH = mesh.NodesTotal(dofRH);
    for (NodeSimple& node : groupNodesRH)
    {
        node.SetValue(0, 1.0);
    }
    Group<NodeSimple> groupNodesWV = mesh.NodesTotal(dofWV);
    for (NodeSimple& node : groupNodesWV)
    {
        node.SetValue(0, 0.18);
    }

    // DOF numbering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DofInfo dofInfo = DofNumbering::Build(mesh.NodesTotal(dofRH), dofRH, constraints);
    dofInfo.Merge(dofWV, DofNumbering::Build(mesh.NodesTotal(dofWV), dofWV, constraints));


    // Create integrand %%%%%%%%%%%%%%%%%%%%%%%%%
    MoistureTransport<1, DCw, DCg, MeC, WVEq> integrandMoistureTransport(dofRH, dofWV, 1., 1., 0.2);
    double delta_t = 0.1;

    auto functionGradientMoistureTransport = std::bind(&MoistureTransport<1, DCw, DCg, MeC, WVEq>::Gradient,
                                                       integrandMoistureTransport, _1, _2, std::ref(delta_t));

    auto functionStiffnessMoistureTransport = std::bind(&MoistureTransport<1, DCw, DCg, MeC, WVEq>::Stiffness,
                                                        integrandMoistureTransport, _1, _2, std::ref(delta_t));

    auto functionDampingMoistureTransport = std::bind(&MoistureTransport<1, DCw, DCg, MeC, WVEq>::Damping,
                                                      integrandMoistureTransport, _1, _2, std::ref(delta_t));

    // Create integration type %%%%%%%%%%%%%%%%%%
    IntegrationTypeTensorProduct<1> integrationType(3, eIntegrationMethod::GAUSS);

    // Create cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    boost::ptr_vector<CellInterface> cellContainer;
    Group<CellInterface> groupCellsMoistureTransport;
    int cellId = 0;
    for (ElementCollection& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, cellId++));
        groupCellsMoistureTransport.Add(cellContainer.back());
    }

    // Assembler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SimpleAssembler assembler(dofInfo);
    GlobalDofVector gradient =
            assembler.BuildVector(groupCellsMoistureTransport, {dofRH, dofWV}, functionGradientMoistureTransport);
    GlobalDofMatrixSparse stiffness =
            assembler.BuildMatrix(groupCellsMoistureTransport, {dofRH, dofWV}, functionStiffnessMoistureTransport);
    GlobalDofMatrixSparse damping =
            assembler.BuildMatrix(groupCellsMoistureTransport, {dofRH, dofWV}, functionDampingMoistureTransport);
}
