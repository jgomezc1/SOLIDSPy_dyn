def respuesta(cells, cell_data, phy_lin):
    """Extracts the nodes located at the physical line
       phy_line

    Parameters
    ----------
        cell : dictionary
            Dictionary created by meshio with cells information.
        cell_data: dictionary
            Dictionary created by meshio with cells data information.
        phy_lin : int
            Physical line to print nodal histories.

    Returns
    -------
        nodes_carga : int
            Array with the nodal data corresponding to the physical
            line phy_line.

    """
    lines = cells["line"]
    phy_line = cell_data["line"]["physical"]
    id_carga = [cont for cont in range(len(phy_line))
                if phy_line[cont] == phy_lin]
    nodes_carga = lines[id_carga]
    nodes_carga = nodes_carga.flatten()
    nodes_carga = list(set(nodes_carga))
    nodes_carga.sort(reverse=False)
    
    return nodes_carga
#
