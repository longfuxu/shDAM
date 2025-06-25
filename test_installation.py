#!/usr/bin/env python3
"""
Installation Test Script for shDAM
Run this script to verify your installation is working correctly.
"""

import sys
import importlib

def test_package(package_name, display_name=None):
    """Test if a package can be imported."""
    if display_name is None:
        display_name = package_name
    
    try:
        importlib.import_module(package_name)
        print(f"✅ {display_name}: OK")
        return True
    except ImportError as e:
        print(f"❌ {display_name}: FAILED - {e}")
        return False

def main():
    """Run installation tests."""
    print("🔍 Testing shDAM Installation")
    print("=" * 40)
    
    # Core packages
    print("\n📦 Core Scientific Packages:")
    core_packages = [
        ("numpy", "NumPy"),
        ("scipy", "SciPy"), 
        ("pandas", "Pandas"),
        ("matplotlib", "Matplotlib"),
    ]
    
    core_success = True
    for pkg, name in core_packages:
        if not test_package(pkg, name):
            core_success = False
    
    # Image processing
    print("\n🖼️ Image Processing:")
    image_packages = [
        ("cv2", "OpenCV"),
        ("tifffile", "TiffFile"),
    ]
    
    image_success = True
    for pkg, name in image_packages:
        if not test_package(pkg, name):
            image_success = False
    
    # Data handling
    print("\n📁 Data File Handling:")
    data_packages = [
        ("nptdms", "NPTDMS"),
        ("openpyxl", "OpenPyXL"),
    ]
    
    data_success = True
    for pkg, name in data_packages:
        if not test_package(pkg, name):
            data_success = False
    
    # Scientific utilities
    print("\n🧮 Scientific Utilities:")
    sci_packages = [
        ("sympy", "SymPy"),
        ("plotly", "Plotly"),
        ("more_itertools", "More Itertools"),
        ("pwlf", "PWLF"),
        ("tabulate", "Tabulate"),
    ]
    
    sci_success = True
    for pkg, name in sci_packages:
        if not test_package(pkg, name):
            sci_success = False
    
    # Optional packages
    print("\n🔬 Optional Packages:")
    optional_packages = [
        ("numba", "Numba (for performance)"),
        ("lumicks.pylake", "Lumicks PyLake"),
    ]
    
    optional_success = 0
    for pkg, name in optional_packages:
        if test_package(pkg, name):
            optional_success += 1
    
    # GUI Framework packages
    print("\n🖥️ GUI Framework Packages:")
    gui_packages = [
        ("customtkinter", "CustomTkinter"),
        ("ttkbootstrap", "TTK Bootstrap"),
        ("seaborn", "Seaborn"),
    ]
    
    gui_success = 0
    for pkg, name in gui_packages:
        if test_package(pkg, name):
            gui_success += 1
    
    # Jupyter
    print("\n📓 Jupyter (Optional):")
    jupyter_packages = [
        ("jupyter", "Jupyter"),
        ("jupyterlab", "JupyterLab"),
        ("ipykernel", "IPython Kernel"),
    ]
    
    jupyter_success = 0
    for pkg, name in jupyter_packages:
        if test_package(pkg, name):
            jupyter_success += 1
    
    # Summary
    print("\n" + "=" * 40)
    print("📊 INSTALLATION SUMMARY")
    print("=" * 40)
    
    if core_success and image_success and data_success and sci_success:
        print("🎉 SUCCESS: All required packages are installed!")
        print("   You can now use the shDAM analysis tools.")
        
        if optional_success > 0:
            print(f"✨ BONUS: {optional_success}/{len(optional_packages)} optional packages available")
        
        if gui_success > 0:
            print(f"🖥️ GUI: {gui_success}/{len(gui_packages)} GUI framework packages available")
            if gui_success == len(gui_packages):
                print("   All GUI tools should work properly!")
        
        if jupyter_success > 0:
            print(f"📓 JUPYTER: {jupyter_success}/{len(jupyter_packages)} Jupyter packages available")
            if jupyter_success == len(jupyter_packages):
                print("   You can run the Jupyter notebooks!")
    else:
        print("⚠️  INCOMPLETE: Some required packages are missing.")
        print("   Please check the installation instructions in README.md")
        
        if sys.platform.startswith('win'):
            print("\n💡 WINDOWS USERS:")
            print("   Consider using conda instead of pip:")
            print("   conda create -n shdam python=3.9 numpy scipy pandas matplotlib opencv jupyter numba")
    
    print("\n" + "=" * 40)
    print(f"Python version: {sys.version}")
    print(f"Platform: {sys.platform}")
    print("=" * 40)

if __name__ == "__main__":
    main() 