# **Find Streets Eulerian Path**  
A Python tool designed to find the shortest path to traverse all street lanes in a city. This is particularly useful for scenarios such as **self-driving car fleet training**, where a human driver must first cover all street lanes before autonomous vehicles can safely operate. By optimizing this path, the tool minimizes **human labor hours** and **fuel or energy costs**.

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
This project is part of a **bachelor's thesis** at the **University of Tartu, Institute of Computer Science**. The primary goal is to develop a tool that models city street networks using graph theory, computes an **Eulerian path** for all connected street lanes, and outputs a path that is optimal for driving through all lanes efficiently. The results are intended for **geospatial visualization in QGIS** and practical applications in optimizing urban logistics.  

While this project was primarily developed for academic purposes, it is publicly available, and contributions are welcome. However, future maintenance may be limited after the thesis defense.

---

## ✨ **Features**
- Builds a graph from **OSM XML** and **GeoJSON** input files representing road networks.
- Balances the graph (if needed) and computes the **Eulerian path** for all connected street segments.
- Exports results to:
  - **Text file**: Eulerian path length, total length of the lanes, and a coefficient (shows how much is the Eulerian path longer than the total length of the lanes).
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

## 🚀 **Usage**

### **Input Files**
The program requires input files generated from the [osm2streets](https://a-b-street.github.io/osm2streets/){:target="_blank"} web app.