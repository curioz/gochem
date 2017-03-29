/*
 * ramachandran.go, part of gochem
 *
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 *
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2.1 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
*/
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chemplot

import (
	"fmt"
	"image/color"
	"math"
	"strings"

	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/vg"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
)

const (
	ErrNilData          = "goChem/ChemPlot: Nil data given "
	ErrInconsistentData = "goChem/ChemPlot: Inconsistent data length "
	ErrTooManyTags      = "goChem/ChemPlot: Maximun number of tagable residues is 4"
	ErrOutOfRange       = "goChem/ChemPlot: Index requested out of range"
)

type Error struct {
	message    string //The error message itself.
	code       string //the name of the QM program giving the problem, or empty string if none
	function   string //the function returning the error.
	additional string //anything else!
	critical   bool
}

func (err Error) Error() string { return fmt.Sprintf("%s  Message: %s", err.function, err.message) }

func (err Error) Code() string { return err.code } //May not be needed

func (err Error) FunctionName() string { return err.function }

func (err Error) Critical() bool { return err.critical }

type RamaSet struct {
	Cprev   int
	N       int
	Ca      int
	C       int
	Npost   int
	MolID   int
	Molname string
}

func basicRamaPlot(title string) (*plot.Plot, error) {
	p, err := plot.New()
	if err != nil {
		return nil, err
	}
	p.Title.Padding = vg.Millimeters(3)
	p.Title.Text = title //"Ramachandran plot"
	p.X.Label.Text = "Phi"
	p.Y.Label.Text = "Psi"
	//Constant axes
	p.X.Min = -180
	p.X.Max = 180
	p.Y.Min = -180
	p.Y.Max = 180
	p.Add(plotter.NewGrid())
	return p, nil

}

// RamaPlotParts produces plots, in png format for the ramachandran data (phi and psi dihedrals)
// contained in data. Data points in tag (maximun 4) are highlighted in the plot.
// the extension must be included in plotname. Returns an error or nil. In RamaPlotParts
// The data is divided in several slices, where each is represented differently in the plot
func RamaPlotParts(data [][][]float64, tag [][]int, title, plotname string) error {
	var err error
	if data == nil {
		return Error{ErrNilData, "", "RamaPlot", "", true}
	}
	// Create a new plot, set its title and
	// axis labels.
	p, err2 := basicRamaPlot(title)
	if err2 != nil {
		return Error{err2.Error(), "", "RamaPlotParts", "", true}
	}
	var tagged int
	for key, val := range data {
		temp := make(plotter.XYs, 1) //len(val))
		//	fmt.Println(key, len(val))
		for k, v := range val {
			temp[0].X = v[0]
			temp[0].Y = v[1]
			// Make a scatter plotter and set its style.
			s, err := plotter.NewScatter(temp) //(pts)
			if err != nil {
				return Error{err.Error(), "", "RamaPlotParts", "", true}

			}
			if tag != nil {
				if len(tag) < len(data) {
					return Error{ErrInconsistentData, "", "RamaPlotParts", "If a non-nil tag slice is provided it must contain an element (which can be nil) for each element in the dihedral slice", true}
				}
				if tag[key] != nil && isInInt(tag[key], k) {
					s.GlyphStyle.Shape, err = getShape(tagged)
					tagged++
				}
			}
			//set the colors
			r, g, b := colors(key, len(data))
			fmt.Println("DATA POINT", key, "color", r, g, b)
			s.GlyphStyle.Color = color.RGBA{R: r, B: b, G: g, A: 255}
			//The tagging procedure is a bit complex.
			p.Add(s)
		}

	}
	filename := fmt.Sprintf("%s.png", plotname)
	//here I  intentionally shadow err.
	if err := p.Save(5, 5, filename); err != nil {
		return Error{err2.Error(), "", "RamaPlotParts", "", true}
	}

	return err
}

//takes hue (0-360), v and s (0-1), returns r,g,b (0-255)
func iHVS2RGB(h, v, s float64) (uint8, uint8, uint8) {
	var i, f, p, q, t float64
	var r, g, b float64
	maxcolor := 255.0
	conversion := maxcolor * v
	if s == 0.0 {
		return uint8(conversion), uint8(conversion), uint8(conversion)
	}
	//conversion:=math.Sqrt(3*math.Pow(maxcolor,2))*v
	h = h / 60
	i = math.Floor(h)
	f = h - i
	p = v * (1 - s)
	q = v * (1 - s*f)
	t = v * (1 - s*(1-f))
	switch int(i) {
	case 0:
		r = v
		g = t
		b = p
	case 1:
		r = q
		g = v
		b = p
	case 2:
		r = p
		g = v
		b = t
	case 3:
		r = p
		g = q
		b = v
	case 4:
		r = t
		g = p
		b = v
	default: //case 5
		r = v
		g = p
		b = q
	}

	r = r * conversion
	g = g * conversion
	b = b * conversion
	return uint8(r), uint8(g), uint8(b)
}

func colors(key, steps int) (r, g, b uint8) {
	norm := 260.0 / float64(steps)
	hp := float64((float64(key) * norm) + 20.0)
	var h float64
	if hp < 55 {
		h = hp - 20.0
	} else {
		h = hp + 20.0
	}
	//	fmt.Println("HUE", h, hp)
	s := 1.0
	v := 1.0
	r, g, b = iHVS2RGB(h, v, s)
	return r, g, b
}

func colorsOld(key, steps int) (r, g, b uint8) {
	norm := (2 * 255.0 / (steps - 1))
	b = uint8(key * norm)
	r = uint8(255) - b
	var critical int
	if norm*(key-1) < 256 && norm*key >= 256 {
		critical = key
	}
	if key*norm > 255 {
		g = uint8(norm * (key - critical))
		b = 255 - g
		r = 0
	}
	//	fmt.Println("crit", critical, norm, steps, key, r, g, b)
	/*	if (key-critical)*norm>255{
			r=uint8(norm*(key-critical))
			g=90
			b=255-r

			}
		}
		fmt.Println(r,g,b, norm, steps)
	*/
	return r, g, b
}

// Produce plots, in png format for the ramachandran data (psi and phi dihedrals)
// contained in data. Data points in tag (maximun 4) are highlighted in the plot.
// the extension must be included in plotname. Returns an error or nil*/
func RamaPlot(data [][]float64, tag []int, title, plotname string) error {
	var err error
	if data == nil {
		return Error{ErrNilData, "", "RamaPlot", "", true}
	}
	// Create a new plot, set its title and
	// axis labels.
	p, err := basicRamaPlot(title)
	if err != nil {
		return Error{err.Error(), "", "RamaPlot", "", true}

	}
	temp := make(plotter.XYs, 1)
	var tagged int //How many residues have been tagged?
	for key, val := range data {
		temp[0].X = val[0]
		temp[0].Y = val[1]
		// Make a scatter plotter and set its style.
		s, err := plotter.NewScatter(temp) //(pts)
		if err != nil {
			return Error{err.Error(), "", "RamaPlot", "", true}
		}
		r, g, b := colors(key, len(data))
		if tag != nil && isInInt(tag, key) {
			s.GlyphStyle.Shape, err = getShape(tagged)
			tagged++
		}
		s.GlyphStyle.Color = color.RGBA{R: r, B: b, G: g, A: 255}
		//		fmt.Println(r,b,g, key, norm, len(data)) //////////////////////////
		// Add the plotter
		p.Add(s)
	}
	// Save the plot to a PNG file.
	filename := fmt.Sprintf("%s.png", plotname)
	//here I  intentionally shadow err.
	if err := p.Save(4, 4, filename); err != nil {
		return Error{err.Error(), "", "RamaPlot", "", true}
	}
	return err
}

func getShape(tagged int) (plot.GlyphDrawer, error) {
	switch tagged {
	case 0:
		return plot.PyramidGlyph{}, nil
	case 1:
		return plot.CircleGlyph{}, nil
	case 2:
		return plot.SquareGlyph{}, nil
	case 3:
		return plot.CrossGlyph{}, nil
	default:
		return plot.RingGlyph{}, Error{ErrTooManyTags, "", "getShape", "", false} // you can still ignore the error and will get just the regular glyph (your residue will not be tagegd)
	}
}

// RamaCalc Obtains the values for the phi and psi dihedrals indicated in []Ramaset, for the
// structure M. The angles are in *degrees*.  It returns a slice of 2-element slices, one for the phi the next for the psi
// dihedral, a and an error or nil.
func RamaCalc(M *v3.Matrix, dihedrals []RamaSet) ([][]float64, error) {
	if M == nil || dihedrals == nil {
		return nil, Error{ErrNilData, "", "RamaCalc", "", true}
	}
	r, _ := M.Dims()
	Rama := make([][]float64, 0, len(dihedrals))
	for _, j := range dihedrals {
		if j.Npost >= r {
			return nil, Error{ErrOutOfRange, "", "RamaCalc", "", true}
		}
		Cprev := M.VecView(j.Cprev)
		N := M.VecView(j.N)
		Ca := M.VecView(j.Ca)
		C := M.VecView(j.C)
		Npost := M.VecView(j.Npost)
		phi := chem.Dihedral(Cprev, N, Ca, C)
		psi := chem.Dihedral(N, Ca, C, Npost)
		temp := []float64{phi * (180 / math.Pi), psi * (180 / math.Pi)}
		Rama = append(Rama, temp)
	}
	return Rama, nil
}

// Filter the set of dihedral angles of a ramachandran plot by residue.(ex. only GLY, everything but GLY)
// The 3 letter code of the residues to be filtered in or out is in filterdata, whether they are filter in
// or out depends on shouldBePresent. It returns the filtered data and a slice containing the indexes in
// the new data of the residues in the old data, when they are included, or -1 when they are not included.
func RamaResidueFilter(dihedrals []RamaSet, filterdata []string, shouldBePresent bool) ([]RamaSet, []int) {
	RetList := make([]RamaSet, 0, 0)
	Index := make([]int, len(dihedrals))
	var added int
	for key, val := range dihedrals {
		isPresent := isInString(filterdata, val.Molname)
		if isPresent == shouldBePresent {
			RetList = append(RetList, val)
			Index[key] = added
			added++
		} else {
			Index[key] = -1
		}
	}
	return RetList, Index
}

// RamaList takes a molecule and returns a slice of RamaSet, which contains the
// indexes for each dihedral to be included in a Ramachandran plot. It gets the dihedral
// indices for all residues in the range resran. if resran has 2 elements defining the
// boundaries. Otherwise, returns dihedral lists for the residues included in
// resran. If resran has 2 elements and the last is -1, RamaList will
// get all the dihedral for residues from resran[0] to the end of the chain.
// It only obtain dihedral lists for residues belonging to a chain included in chains.
func RamaList(M chem.Atomer, chains string, resran []int) ([]RamaSet, error) {
	RamaList := make([]RamaSet, 0, 0)
	if len(resran) == 2 {
		if resran[1] == -1 {
			resran[1] = 999999999 //should work!
		}
	}
	if M == nil {
		return nil, Error{ErrNilData, "", "RamaList", "nil molecule", true}
	}
	C := -1
	N := -1
	Ca := -1
	Cprev := -1
	Npost := -1
	chainprev := "NOTAVALIDCHAIN" //any non-valid chain name
	for num := 0; num < M.Len(); num++ {
		at := M.Atom(num)
		//First get the indexes we need. Change: If you give RamaList an empty string for "chains", it will 
		//include all chains in the chem.Atomer.
		if strings.Contains(chains, string(at.Chain)) || at.Chain == " " || chains=="" {
			if at.Chain != chainprev {
				chainprev = at.Chain
				C = -1
				N = -1
				Ca = -1
				Cprev = -1
				Npost = -1
			}
			if at.Name == "C" && Cprev == -1 {
				Cprev = num
			}
			if at.Name == "N" && Cprev != -1 && N == -1 && at.MolID > M.Atom(Cprev).MolID {
				N = num
			}
			if at.Name == "C" && Cprev != -1 && at.MolID > M.Atom(Cprev).MolID {
				C = num
			}
			if at.Name == "CA" && Cprev != -1 && at.MolID > M.Atom(Cprev).MolID {
				Ca = num
			}
			if at.Name == "N" && Ca != -1 && at.MolID > M.Atom(Ca).MolID {
				Npost = num
			}
			//when we have them all, we save
			if Cprev != -1 && Ca != -1 && N != -1 && C != -1 && Npost != -1 {
				//We check that the residue ids are what they are supposed to be
				r1 := M.Atom(Cprev).MolID
				r2 := M.Atom(N).MolID
				r2a := M.Atom(Ca).MolID
				r2b := M.Atom(C).MolID
				r3 := M.Atom(Npost).MolID
				if (len(resran) == 2 && (r2 >= resran[0] && r2 <= resran[1])) || isInInt(resran, r2) {
					if r1 != r2-1 || r2 != r2a || r2a != r2b || r2b != r3-1 {
						return nil, Error{fmt.Sprintf("Incorrect backbone Cprev: %d N-1: %d CA: %d C: %d Npost-1: %d", r1, r2-1, r2a, r2b, r3-1), "", "RamaList", "", true}
					}
					temp := RamaSet{Cprev, N, Ca, C, Npost, r2, M.Atom(Ca).Molname}
					RamaList = append(RamaList, temp)
				}
				N = Npost
				Ca = -1
				Cprev = C
				C = -1
				Npost = -1
			}
		}
	}
	//	fmt.Println("Rama",Rama, "failed", failed)
	return RamaList, nil
}
