# Start with a Julia base image
FROM julia:1.11.3

# Set up the Julia environment
WORKDIR /app

# Install just the computational packages. BuildCGAL downloads its prebuilt
# CGAL tools here, so no compiler, Boost or CGAL build is needed in the image.
RUN julia -e 'using Pkg; \
    Pkg.add(url="https://github.com/Timmmdavis/FischerNewton.git"); \
    Pkg.add(url="https://github.com/Timmmdavis/CutAndDisplaceJulia.git"); \
    Pkg.add(url="https://github.com/Timmmdavis/BuildCGAL.git"); \
    Pkg.add("Revise"); \
    Pkg.precompile()'

# Check the CGAL tools are installed and run on sample meshes
RUN julia -e 'using Pkg; Pkg.test("BuildCGAL")'

# Create a test script that demonstrates the functionality
RUN echo 'println("=== CutAndDisplaceJulia Test Environment ===");\
         using FischerNewton;\
         using CutAndDisplaceJulia;\
         using BuildCGAL;\
         using Revise;\
         println("\\nAll core packages loaded successfully!");\
         println("\\nNOTE: CutAndDisplaceJuliaPlots is not available in this container.");\
         println("      For visualization, save your results to the shared volume and");\
         println("      view them with external tools like Meshlab.\\n");\
         println("Ready to use CutAndDisplaceJulia!")' > /app/startup.jl

# Create a helper script to run the dyke propagation test
RUN echo 'function run_dyke_propagation(output_dir="/data")\
    cd(output_dir);\
    ModuleDir = dirname(dirname(pathof(CutAndDisplaceJulia)));\
    TestScript = joinpath(ModuleDir, "test", "test_PropagationTest.jl");\
    include(TestScript);\
    # Configuration\
    HSFlag = 0;  # Constant\
    printstyled("Hs off \\n", color=:cyan);\
    # Elastic constants\
    G = ShearModulus(2.0e9);\
    ν = PoissonsRatio(0.25);\
    (K, E, λ, ν, G) = CutAndDisplaceJulia.ElasticConstantsCheck(G, ν);\
    # Material properties\
    ρrock = 2900;\
    ρfluid = 2600;\
    g = 9.81;\
    Δρ = ((ρrock-ρfluid)*g);\
    Kc = 5e7;  # 50 MPa √m\
    # Crack volume calculation\
    CrackVolume = real(((1-ν)/(16*G))*((9*π^4*Kc^8)/(complex(Δρ)^5))^(1/3));\
    NoTris = 300;\
    # Run simulation\
    (PropFlag, maxX, minX, maxY, minY, maxZ, minZ) = testProp(HSFlag, ν, G, Δρ, Kc, CrackVolume, NoTris);\
    println("\\nDyke propagation test completed. Results saved to: $output_dir");\
end' > /app/dyke_test.jl

# Create README with usage instructions in the container
RUN echo '# CutAndDisplaceJulia Docker Container\n\
\n\
This container includes the computational components of CutAndDisplaceJulia with prebuilt CGAL tools (via BuildCGAL).\n\
\n\
## Usage\n\
\n\
### Basic Julia REPL with packages loaded:\n\
```\n\
docker run -it --rm -v /path/to/local/folder:/data cutanddisplace-julia\n\
```\n\
\n\
### Running the dyke propagation test:\n\
```julia\n\
include("/app/dyke_test.jl")\n\
run_dyke_propagation("/data")\n\
```\n\
\n\
### Viewing results:\n\
Results are saved to your mounted volume. Use external tools like Meshlab to view mesh files.\n\
\n\
### Note about visualization:\n\
CutAndDisplaceJuliaPlots is not available in this container due to its dependency on\n\
OpenGL, which is difficult to run in a headless Docker environment. For visualization,\n\
save your results to the mounted volume and use external tools.\n\
' > /app/README.md

ENTRYPOINT ["julia", "-L", "/app/startup.jl"]