/*
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <http://unlicense.org>
 */

#include "geant4.hh"
/* For the run and UI managers */
#include "G4RunManager.hh"
/* For the detector construction */
#include "G4Ellipsoid.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4Turtle.hh"
#include "G4PVPlacement.hh"
/* For the Physics */
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronicProcessStore.hh"
/* For the primary generator */
#include "G4Geantino.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ParticleGun.hh"
/* For the random engine */
#include "Randomize.hh"
/* C interfaces */
extern "C" {
#include "downsampler.h"
#include "geant4.h"
#include "geometry.h"
#include "turtle.h"
}
/* Unix tricks */
#include <fcntl.h>
#include <unistd.h>


#ifndef INSTALL_PREFIX
#define INSTALL_PREFIX "."
#endif


static G4RunManager * gRun = NULL;
static G4int gEvents = 1;

static struct {
        G4ParticleDefinition * particle;
        G4double energy;
        G4ThreeVector position;
        G4ThreeVector direction;
} gPrimary;

static struct {
        G4String area;
        int period;
        enum geant4_mode mode;
} gTopography;

static struct {
        G4int steps;
        G4double depth;
        G4double time;
} gStats;


static struct turtle_map * gMap, * gGeoid;
static const struct turtle_projection * gProjection = NULL;


static void toggle_stdout()
{
        static int bak, muted = 0;

        if (muted) {
                fflush(stdout);
                dup2(bak, 1);
                close(bak);
                muted = 0;
        } else {
                fflush(stdout);
                bak = dup(1);
                const int tmp = open("/dev/null", O_WRONLY);
                dup2(tmp, 1);
                close(tmp);
                muted = 1;
        }
}


static G4double GetZMax()
{
        if (gTopography.area == "tianshan")
                return GEOMETRY_ZMAX_TIANSHAN * CLHEP::m;
        else
                return GEOMETRY_ZMAX_PDD * CLHEP::m;
}


static G4ThreeVector GetNode(int ix, int iy)
{
        double x, y, z;
        turtle_map_node(gMap, ix, iy, &x, &y, &z);

        if (gProjection != NULL) {
                turtle_projection_unproject(
                    gProjection, x, y, &y, &x);
        }

        double undulation;
        turtle_map_elevation(gGeoid, x, y, &undulation, NULL);
        z += undulation;

        double ecef[3];
        turtle_ecef_from_geodetic(y, x, z, ecef);
        return G4ThreeVector(
            ecef[0] * CLHEP::m, ecef[1] * CLHEP::m, ecef[2] * CLHEP::m);
}


static G4ThreeVector GetBottomNode(
    int ix, int iy, const G4ThreeVector & vertical)
{
        G4ThreeVector n = GetNode(ix, iy);

        double altitude;
        double ecef[3] = { n.x() / CLHEP::m, n.y() / CLHEP::m,
                           n.z() / CLHEP::m };
        turtle_ecef_to_geodetic(ecef, NULL, NULL, &altitude);

        n -= (altitude + 100.0) * CLHEP::m * vertical;
        return n;
}


static void WriteFacet(FILE * fid, const G4ThreeVector& n0,
   const G4ThreeVector& n1, const G4ThreeVector& n2)
{
        float r[13] = { 0., 0., 0.,
                       (float)n0.x(), (float)n0.y(), (float)n0.z(),
                       (float)n1.x(), (float)n1.y(), (float)n1.z(),
                       (float)n2.x(), (float)n2.y(), (float)n2.z(), 0 };

        fwrite(r, sizeof(char), 50, fid);
}


static void AddQuad(G4TessellatedSolid * geometry, FILE * fid,
     const G4ThreeVector& n00, const G4ThreeVector& n10,
     const G4ThreeVector& n01, const G4ThreeVector& n11)
{
        G4TriangularFacet * facet;
        facet = new G4TriangularFacet(n00, n10, n11, ABSOLUTE);
        geometry->AddFacet((G4VFacet *)facet);
        facet = new G4TriangularFacet(n11, n01, n00, ABSOLUTE);
        geometry->AddFacet((G4VFacet *)facet);

        if (fid != NULL) {
                WriteFacet(fid, n00, n10, n11);
                WriteFacet(fid, n11, n01, n00);
        }
}


static G4VPhysicalVolume * BuildTessellation()
{
        /* Load the map data */
        G4String path = G4String(INSTALL_PREFIX) + "/share/topography/" +
                        gTopography.area + ".png";
        downsampler_map_load(&gMap, path.c_str(), gTopography.period);
        gProjection = turtle_map_projection(gMap);
        turtle_map_load(&gGeoid, GEOMETRY_GEOID);

        struct turtle_map_info info;
        turtle_map_meta(gMap, &info, NULL);
        const int nx = info.nx, ny = info.ny;

        /* Initialise the solid */
        G4TessellatedSolid * geometry =
            new G4TessellatedSolid("Topography");

        FILE * fid = NULL;
        long int pos = 0;
        if (gTopography.mode == GEANT4_MODE_STL) {
                char buffer[80];
                const int nh = sizeof(buffer) / sizeof(*buffer);
                snprintf(buffer, nh, "%s-%d.stl",
                    gTopography.area.c_str(), gTopography.period);

                fid = fopen(buffer, "wb+");
                if (fid == NULL) {
                        fprintf(stderr, "could not create STL file\n");
                        exit(EXIT_FAILURE);
                }

                memset(buffer, 0x0, nh);
                snprintf(buffer, nh, "%s %d",
                    gTopography.area.c_str(), gTopography.period);
                fwrite(buffer, sizeof *buffer, nh, fid);
                pos = ftell(fid);
                uint32_t tmp = 0;
                fwrite(&tmp, sizeof tmp, 1, fid);
        }

        /* Build the facets */
        std::cout << "# building facets ...";
        struct tictac tictac;
        tictac_initialise(&tictac);
        tictac.start(&tictac);

        int n_iter = (nx - 1) * (ny - 1);
        int period = n_iter / 100;
        if (period <= 0) period = 1;

        int ix, iy, it = 0;
        for (iy = 0; iy < ny - 1; iy++) for (ix = 0; ix < nx - 1; ix++, it++) {
                G4ThreeVector n00 = GetNode(ix, iy);
                G4ThreeVector n10 = GetNode(ix + 1, iy);
                G4ThreeVector n01 = GetNode(ix, iy + 1);
                G4ThreeVector n11 = GetNode(ix + 1, iy + 1);

                AddQuad(geometry, fid, n00, n10, n01, n11);

                if ((it + 1) % period == 0) {
                        std::cout << "\r# building facets ... "
                                  << (it + 1) / period << " %";
                        fflush(stdout);
                }
        }

        /* Build the bottom facets */
        G4ThreeVector origin = GetNode(nx / 2, ny / 2);
        double latitude0, longitude0;
        double ecef[3] = { origin.x() / CLHEP::m, origin.y() / CLHEP::m,
                           origin.z() / CLHEP::m };
        turtle_ecef_to_geodetic(ecef, &latitude0, &longitude0, NULL);
        turtle_ecef_from_horizontal(latitude0, longitude0, 0, 90, ecef);
        G4ThreeVector vertical(ecef[0], ecef[1], ecef[2]);

        G4ThreeVector b00 = GetBottomNode(0, 0, vertical);
        G4ThreeVector b10 = GetBottomNode(nx - 1, 0, vertical);
        G4ThreeVector b01 = GetBottomNode(0, ny - 1, vertical);
        G4ThreeVector b11 = GetBottomNode(nx - 1, ny - 1, vertical);

        AddQuad(geometry, fid, b11, b10, b01, b00);
        n_iter++;

        /* Create the border facets */
        for (iy = 0; iy < ny; iy += ny - 1) {
                G4ThreeVector & b0 = iy ? b01 : b00;
                G4ThreeVector & b1 = iy ? b11 : b10;

                for (ix = 0; ix < nx - 1; ix++, n_iter++) {
                        const double h0 = ix / float(nx - 1);
                        const double h1 = (ix + 1) / float(nx - 1);

                        G4ThreeVector n00 = GetNode(ix, iy);
                        G4ThreeVector n10 = GetNode(ix + 1, iy);
                        G4ThreeVector n01 = b0 * (1. - h0) + b1 * h0;
                        G4ThreeVector n11 = b0 * (1. - h1) + b1 * h1;

                        if (iy) AddQuad(geometry, fid, n00, n10, n01, n11);
                        else    AddQuad(geometry, fid, n11, n10, n01, n00);
                }
        }

        for (ix = 0; ix < nx; ix += nx - 1) {
                G4ThreeVector & b0 = ix ? b10 : b00;
                G4ThreeVector & b1 = ix ? b11 : b01;

                for (iy = 0; iy < ny - 1; iy++, n_iter++) {
                        const double h0 = iy / float(ny - 1);
                        const double h1 = (iy + 1) / float(ny - 1);

                        G4ThreeVector n00 = GetNode(ix, iy);
                        G4ThreeVector n10 = GetNode(ix, iy + 1);
                        G4ThreeVector n01 = b0 * (1. - h0) + b1 * h0;
                        G4ThreeVector n11 = b0 * (1. - h1) + b1 * h1;

                        if (ix) AddQuad(geometry, fid, n11, n10, n01, n00);
                        else    AddQuad(geometry, fid, n00, n10, n01, n11);
                }
        }

        double dt = tictac.stop(&tictac);
        std::cout << "\r# " <<  2 * n_iter << " facets built in "
                  << dt << " s" << std::endl;

        /* Free the maps */
        turtle_map_destroy(&gMap);
        turtle_map_destroy(&gGeoid);
        gProjection = NULL;

        if (fid != NULL) {
                fseek(fid, pos, SEEK_SET);
                uint32_t tmp = 2 * n_iter;
                fwrite(&tmp, sizeof tmp, 1, fid);
                fclose(fid);
        }

        /* Finalise the geometry */
        toggle_stdout();
        tictac.start(&tictac);
        geometry->SetSolidClosed(true);
        dt = tictac.stop(&tictac);
        toggle_stdout();
        std::cout << "\r# topography closed in " << dt << " s" << std::endl;

        /* Create the scene */
        const double hmax = GetZMax();
        const double a = 6378137.0 * CLHEP::m + hmax;
        const double b = 6356752.3142 * CLHEP::m + hmax;
        G4Ellipsoid * worldGeometry = new G4Ellipsoid("World", a, a, b);

        G4Turtle * turtle = G4Turtle::GetInstance();
        G4Material * air = turtle->GetAirMaterial();
        G4Material * rock = turtle->GetRockMaterial();

        G4LogicalVolume * worldLogical =
            new G4LogicalVolume(worldGeometry, air, "logical", 0, 0, 0);
        G4VPhysicalVolume * worldPhysical = new G4PVPlacement(0,
            G4ThreeVector(0., 0., 0.), worldLogical, "World", 0, false, 0, 0);

        G4LogicalVolume * topoLogical =
            new G4LogicalVolume(geometry, rock, "logical", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), topoLogical,
            "Topography", worldLogical, false, 0, 0);

        return worldPhysical;
}


/* Build the TURTLE geometry */
static G4VPhysicalVolume * BuildTurtle()
{
        /* Generate the downsampling, if needed */
        G4String path = G4String(INSTALL_PREFIX) + "/share/topography/" +
                        gTopography.area + ".png";
        if (gTopography.period > 1) {
                struct turtle_map * map;
                downsampler_map_load(&map, path.c_str(), gTopography.period);

                path = G4String("/tmp/geant4.topography.") + gTopography.area +
                       ".png";
                turtle_map_dump(map, path.c_str());
                turtle_map_destroy(&map);
        }

        /* Configure the topography */
        G4Turtle * turtle = G4Turtle::GetInstance();
        turtle->SetTopographyData("", path,
            INSTALL_PREFIX "/share/topography/egm96.png");
        turtle->SetBottomLevel(0. * CLHEP::m);
        if (gTopography.area == "tianshan")
                turtle->SetTopLevel(GEOMETRY_ZMAX_TIANSHAN * CLHEP::m);
        else
                turtle->SetTopLevel(GEOMETRY_ZMAX_PDD * CLHEP::m);

        /* Return the topography as world volume */
        return turtle->GetPhysicalVolume();
}


/* Wrap the geometry builder */
G4VPhysicalVolume * DetectorConstruction::Construct()
{
        if (gTopography.mode == GEANT4_MODE_G4TURTLE) {
                return BuildTurtle();
        } else {
                return BuildTessellation();
        }
}


PhysicsList::PhysicsList(G4double cut) : G4VModularPhysicsList()
{
        G4HadronicProcessStore::Instance()->SetVerbose(0);

        /* Set the cut value */
        defaultCutValue = cut;

        /* Instanciate the default Physics constructors */
        fDecayPhysics = new G4DecayPhysics();
        fEmPhysics = new G4EmStandardPhysics();
        fExtraPhysics = new G4EmExtraPhysics();
}


PhysicsList::~PhysicsList()
{
        delete fDecayPhysics;
        delete fEmPhysics;
        delete fExtraPhysics;
}


void PhysicsList::ConstructParticle()
{
        fDecayPhysics->ConstructParticle();
}


void PhysicsList::ConstructProcess()
{
        AddTransportation();
        fDecayPhysics->ConstructProcess();
        fEmPhysics->ConstructProcess();
        fExtraPhysics->ConstructProcess();
}


void PhysicsList::SetCuts()
{
        G4int temp = GetVerboseLevel();
        SetVerboseLevel(0);
        SetCutsWithDefault();
        SetVerboseLevel(temp);
}


PrimaryGenerator::PrimaryGenerator() : G4VUserPrimaryGeneratorAction()
{
        fParticleGun = new G4ParticleGun(1);
}


PrimaryGenerator::~PrimaryGenerator()
{
        delete fParticleGun;
        fParticleGun = NULL;
}


void PrimaryGenerator::GeneratePrimaries(G4Event * event)
{
        fParticleGun->SetParticleDefinition(gPrimary.particle);
        fParticleGun->SetParticleEnergy(gPrimary.energy);
        fParticleGun->SetParticlePosition(gPrimary.position);
        fParticleGun->SetParticleMomentumDirection(gPrimary.direction);
        fParticleGun->GeneratePrimaryVertex(event);
}


G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(
    const G4Track * track)
{
        G4ClassificationOfNewTrack classification = fWaiting;
        G4ParticleDefinition * particleType = track->GetDefinition();

        if ((particleType == G4Geantino::Definition()) ||
            (particleType == G4MuonPlus::Definition()) ||
            (particleType == G4MuonMinus::Definition()))
                classification = fUrgent;
        else
                classification = fKill;

        return classification;
}


RunAction::RunAction() : G4UserRunAction()
{
        tictac_initialise(&fTictac);
}


void RunAction::BeginOfRunAction(const G4Run *)
{
        fTictac.start(&fTictac);
}


void RunAction::EndOfRunAction(const G4Run *)
{
        gStats.time = fTictac.stop(&fTictac);
}


void SteppingAction::UserSteppingAction(const G4Step * step)
{
        static const G4Material * rock =
            G4Turtle::GetInstance()->GetRockMaterial();

        gStats.steps++;

        if ((step->GetTrack()->GetTrackID() > 1) ||
            (step->GetPreStepPoint()->GetMaterial() != rock))
                return;

        gStats.depth += step->GetStepLength();
}


static void geant4_clear(void)
{
        delete gRun;
}


static int geant4_run(enum geant4_particle particle, double energy,
    double position[3], double direction[3], double * depth, double * time)
{
        if (particle == GEANT4_PARTICLE_GEANTINO)
                gPrimary.particle = G4Geantino::Definition();
        else if (particle == GEANT4_PARTICLE_MUON_MINUS)
                gPrimary.particle = G4MuonPlus::Definition();
        else if (particle == GEANT4_PARTICLE_MUON_PLUS)
                gPrimary.particle = G4MuonPlus::Definition();
        else {
                fprintf(stderr, "error: unsuported geant4 particle\n");
                exit(EXIT_FAILURE);
        }

        gPrimary.energy = energy * CLHEP::GeV;
        gPrimary.position = G4ThreeVector(position[0] * CLHEP::m,
            position[1] * CLHEP::m, position[2] * CLHEP::m);
        gPrimary.direction = G4ThreeVector(direction[0] * CLHEP::m,
            direction[1] * CLHEP::m, direction[2] * CLHEP::m);

        memset(&gStats, 0x0, sizeof gStats);
        gRun->BeamOn(gEvents);

        *depth = gStats.depth / (gEvents * CLHEP::m);
        *time = gStats.time / gEvents;

        return gStats.steps / gEvents;
}


extern "C" void geant4_initialise(struct geant4 * geant4,
    const char * topography, int period, enum geant4_mode mode,
    int kill_secondaries, double cut)
{
        gTopography.area = topography;
        gTopography.period = period;
        gTopography.mode = mode;


        /* Set the random seed */
        unsigned long seed = 0;
        FILE * fid = fopen("/dev/urandom", "rb");
        fread(&seed, sizeof(long), 1, fid);
        fclose(fid);
        CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
        G4Random::setTheSeed(seed);

        /* Instanciate the run manager */
        toggle_stdout();
        gRun = new G4RunManager;
        toggle_stdout();

        gRun->SetUserInitialization(new DetectorConstruction);
        gRun->SetUserInitialization(new PhysicsList(cut * CLHEP::m));
        gRun->SetUserAction(new PrimaryGenerator);
        if (kill_secondaries)
                gRun->SetUserAction(new StackingAction);
        gRun->SetUserAction(new RunAction);
        gRun->SetUserAction(new SteppingAction);
        gRun->Initialize();

        geant4->clear = &geant4_clear;
        geant4->run = &geant4_run;
}


extern "C" void geant4_configure(
    struct geant4 *, int events, int verbosity)
{
        if (events > 0) {
                gEvents = events;
        }

        if (verbosity >= 0) {
                G4EventManager * event = G4EventManager::GetEventManager();
                G4TrackingManager * tracking = event->GetTrackingManager();
                tracking->SetVerboseLevel(verbosity);
        }
}
