/*
 * corrections.go, part of gochem.
 *
 *
 * Copyright 2013 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 *
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package qm

import (
	"bufio"
	"fmt"
	"github.com/rmera/gochem"
	"os"
	"os/exec"
	"strconv"
	"strings"
)



//Correction obtains corrections to the energy (in the future probably also to the gradients
//From an extenal program. The executable for the program to be used used (for instance,
//dftd3 from Stefan Grimme) must be in the path
type Correction struct {
	program    string
	damping    string
	functional string
	corrVal    float64
}

//Returns a new Correction structure with the default choices, and with the kind of correction given.
func NewCorrection(kind string) (*Correction, error) {
	C := new(Correction)
	err := C.SetCorrection(kind)
	if err != nil {
		return nil, err
	}
	C.SetDefaults()
	return C, nil
}

//Set the type of correction to be used
func (C *Correction) SetCorrection(kind string) error {
	kind = strings.ToLower(kind)
	switch kind {
	case "d3":
		C.program = "dftd3"

	default:
		return fmt.Errorf("%s correction is not supported", kind)
	}
	return nil
}

func (C *Correction) SetDefaults() {
	C.damping = "-zero"
}

//Sets the Correction damping. "zero" is zero damping, anything else is BJ damping
func (C *Correction) SetDamping(d string) {
	nd := strings.ToLower(d)
	switch nd {
	case "zero":
		C.damping = "-zero"
	default:
		C.damping = "-bj"
	}
}

//Sets the functional to be used in the correction
func (C *Correction) SetFunctional(f string) {
	C.functional = tMMethods[f]
}

//Returns the type of correction currently being performed
func (C *Correction) Correction() string {
	switch C.program {
	case "dftd3":
		return "D3"
	default:
		return ""
	}
}

//Returns the correction to the energy and an error or nil.
func (C *Correction) CorrectE(coords *chem.VecMatrix, ref chem.Atomer) (float64, error) {
	err := chem.XYZWrite("temp.xyz", coords, ref)
	if err != nil {
		return 0, err
	}
	var corr float64
	switch C.program {
	case "dftd3":
		corr, err = C.runD3()

	}
	return corr, err
}


//Runs the dftd3 program from Stefan Grimme, which must be in the path
//and retrieves the correction in kcal/mol.
func (C *Correction) runD3() (float64, error) {
	var energy float64
	out, err := os.Create(fmt.Sprintf("%s.out", C.program))
	if err != nil {
		return 0, err
	}
	command := exec.Command(C.program, "temp.xyz", "-func", C.functional, C.damping, "-anal")
	command.Stdout = out
	err = command.Run()
	if err != nil {
		out.Close()
		return 0, err
	}
	out.Close()
	dispT, err := os.Open(fmt.Sprintf("%s.out", C.program))
	if err != nil {
		return 0, err
	}
	defer dispT.Close()
	disp := bufio.NewReader(dispT)
	for {
		line, err := disp.ReadString('\n')
		if err != nil {
			break
		}
		if strings.Contains(line, "Edisp /kcal,au") {
			energy, err = strconv.ParseFloat(line[17:28], 64)
			break
		}
	}
	return energy, err
}
