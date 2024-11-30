# **Find Streets Eulerian Path**  
A Python tool to generate Eulerian paths for city street networks using OSM (OpenStreetMap) and GeoJSON files. The program balances a graph representing city streets, finds the optimal Eulerian path, and exports the results for further visualization in QGIS.

---

## 📋 **Table of Contents**
- [About the Project](#about-the-project)
- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Installing Dependencies](#installing-dependencies)
- [Usage](#usage)
  - [Input Files](#input-files)
  - [Running the Program](#running-the-program)
- [Output](#output)
- [Contributing](#contributing)
- [License](#license)

---

## 🏛️ **About the Project**
This project is part of a **bachelor's thesis** at the **University of Tartu, Institute of Computer Science**. The goal is to create a tool that models city street networks using graph theory and computes an **Eulerian path** through all street segments, making the results usable for geospatial visualization in QGIS.  

This project was primarily developed for academic purposes, no contributions will be accepted until the thesis has been defended.

---

## ✨ **Features**
- Builds a graph from **OSM XML** and **GeoJSON** input files representing road networks.
- Balances the graph and computes the **Eulerian path** for all connected street segments.
- Exports results to:
  - **Text file**: Eulerian path length, total street length, and a coefficient comparing both.
  - **GeoPackage (GPKG)**: Nodes and edges for easy visualization and animation in **QGIS**.

---

## ⚙️ **Installation**

### **Prerequisites**
Ensure you have Python 3.13 or higher installed. You’ll also need `pip` for package management.

### **Installing Dependencies**
1. Clone this repository:
    ```bash
   git clone https://github.com/raunoraa/find-streets-eulerian-path.git
   cd find-streets-eulerian-path
    ```
2. Install the required Python packages:
    ```bash
   pip install -r requirements.txt
    ```

---