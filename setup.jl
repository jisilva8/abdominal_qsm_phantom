# Setup script for Abdominal QSM Phantom toolbox
# Activates project and installs dependenciesusing Pkg

using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

println("✅ Setup completed: Abdominal QSM Phantom toolbox is ready to use.")