// Package bn256 implements a particular bilinear group.
//
// Bilinear groups are the basis of many of the new cryptographic protocols that
// have been proposed over the past decade. They consist of a triplet of groups
// (G₁, G₂ and GT) such that there exists a function e(g₁ˣ,g₂ʸ)=gTˣʸ (where gₓ
// is a generator of the respective group). That function is called a pairing
// function.
//
// This package specifically implements the Optimal Ate pairing over a 256-bit
// Barreto-Naehrig curve as described in
// http://cryptojedi.org/papers/dclxvi-20100714.pdf. Its output is compatible
// with the implementation described in that paper.
//
// This package previously claimed to operate at a 128-bit security level.
// However, recent improvements in attacks mean that is no longer true. See
// https://moderncrypto.org/mail-archive/curves/2016/000740.html.
package bn256

import (
	"crypto/rand"
	"crypto/sha256"
	"errors"
	"io"
	"math/big"
	"fmt"
)

func randomK(r io.Reader) (k *big.Int, err error) {
	for {
		k, err = rand.Int(r, Order)
		if k.Sign() > 0 || err != nil {
			return
		}
	}

	return
}

// G1 is an abstract cyclic group. The zero value is suitable for use as the
// output of an operation, but cannot be used as an input.
type G1 struct {
	p *curvePoint
}

// RandomG1 returns x and g₁ˣ where x is a random, non-zero number read from r.
func RandomG1(r io.Reader) (*big.Int, *G1, error) {
	k, err := randomK(r)
	if err != nil {
		return nil, nil, err
	}

	return k, new(G1).ScalarBaseMult(k), nil
}

// HashG1 hashes string m to an element in group G1 using
// try and increment method.
func HashG1(m string) (*G1, error) {
	h := sha256.Sum256([]byte(m))
	hashNum := new(big.Int)
	for {
		hashNum.SetBytes(h[:])
		if hashNum.Cmp(p) == -1 {
			break
		}
		h = sha256.Sum256(h[:])
	}

	x, x2, x3, rhs, y := &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}
	x = x.SetInt(hashNum)
	three := newGFp(3)
	var err error
	for {
		//let's check if there exists a point (X, Y) for some Y on EC -
		// that means X^3 + 3 needs to be a quadratic residue

		gfpMul(x2, x, x)
		gfpMul(x3, x2, x)
		gfpAdd(rhs, x3, three)

		y, err = y.Sqrt(rhs) // TODO: what about -y
		if err == nil {      // alternatively, if Y is not needed, big.Jacobi(rhs, P) can be used to check if rhs is quadratic residue
			// BN curve has cofactor 1 (all points of the curve form a group where we are operating),
			// so X (now that we know rhs is QR) is an X-coordinate of some point in a cyclic group
			point := &curvePoint{
				x: *x,
				y: *y,
				z: *newGFp(1),
				t: *newGFp(1),
			}
			return &G1{point}, nil
		}
		gfpAdd(x, x, newGFp(1))
	}
}

func (g *G1) String() string {
	return "bn256.G1" + g.p.String()
}

// ScalarBaseMult sets e to g*k where g is the generator of the group and then
// returns e.
func (e *G1) ScalarBaseMult(k *big.Int) *G1 {
	if e.p == nil {
		e.p = &curvePoint{}
	}
	e.p.Mul(curveGen, k)
	return e
}

// ScalarMult sets e to a*k and then returns e.
func (e *G1) ScalarMult(a *G1, k *big.Int) *G1 {
	if e.p == nil {
		e.p = &curvePoint{}
	}
	e.p.Mul(a.p, k)
	return e
}

// Add sets e to a+b and then returns e.
func (e *G1) Add(a, b *G1) *G1 {
	if e.p == nil {
		e.p = &curvePoint{}
	}
	e.p.Add(a.p, b.p)
	return e
}

// Neg sets e to -a and then returns e.
func (e *G1) Neg(a *G1) *G1 {
	if e.p == nil {
		e.p = &curvePoint{}
	}
	e.p.Neg(a.p)
	return e
}

// Set sets e to a and then returns e.
func (e *G1) Set(a *G1) *G1 {
	if e.p == nil {
		e.p = &curvePoint{}
	}
	e.p.Set(a.p)
	return e
}

// Marshal converts e to a byte slice.
func (e *G1) Marshal() []byte {
	// Each value is a 256-bit number.
	const numBytes = 256 / 8

	if e.p == nil {
		e.p = &curvePoint{}
	}

	e.p.MakeAffine()
	ret := make([]byte, numBytes*2)
	if e.p.IsInfinity() {
		return ret
	}
	temp := &gfP{}

	montDecode(temp, &e.p.x)
	temp.Marshal(ret)
	montDecode(temp, &e.p.y)
	temp.Marshal(ret[numBytes:])

	return ret
}

// Unmarshal sets e to the result of converting the output of Marshal back into
// a group element and then returns e.
func (e *G1) Unmarshal(m []byte) ([]byte, error) {
	// Each value is a 256-bit number.
	const numBytes = 256 / 8

	if len(m) < 2*numBytes {
		return nil, errors.New("bn256: not enough data")
	}

	if e.p == nil {
		e.p = &curvePoint{}
	} else {
		e.p.x, e.p.y = gfP{0}, gfP{0}
	}

	e.p.x.Unmarshal(m)
	e.p.y.Unmarshal(m[numBytes:])
	montEncode(&e.p.x, &e.p.x)
	montEncode(&e.p.y, &e.p.y)

	zero := gfP{0}
	if e.p.x == zero && e.p.y == zero {
		// This is the point at infinity.
		e.p.y = *newGFp(1)
		e.p.z = gfP{0}
		e.p.t = gfP{0}
	} else {
		e.p.z = *newGFp(1)
		e.p.t = *newGFp(1)

		if !e.p.IsOnCurve() {
			return nil, errors.New("bn256: malformed point")
		}
	}

	return m[2*numBytes:], nil
}

// G2 is an abstract cyclic group. The zero value is suitable for use as the
// output of an operation, but cannot be used as an input.
type G2 struct {
	p *twistPoint
}

// RandomG2 returns x and g₂ˣ where x is a random, non-zero number read from r.
func RandomG2(r io.Reader) (*big.Int, *G2, error) {
	k, err := randomK(r)
	if err != nil {
		return nil, nil, err
	}

	return k, new(G2).ScalarBaseMult(k), nil
}

// HashG2 hashes string m to an element in group G2. It uses:
// Fuentes-Castaneda, Laura, Edward Knapp, and Francisco Rodríguez-Henríquez. "Faster hashing to G_2."
// International Workshop on Selected Areas in Cryptography. Springer, Berlin, Heidelberg, 2011.
func HashG2(m string) (*G2, error) {
	h := sha256.Sum256([]byte(m))
	hashNum := new(big.Int)
	for {
		hashNum.SetBytes(h[:])
		if hashNum.Cmp(p) == -1 {
			break
		}
		h = sha256.Sum256(h[:])
	}
	v := &gfP{}
	v.SetInt(hashNum)

	// gfp2 is (x1, y1) where x1*i + y1
	x, xxx, rhs, y := &gfP2{}, &gfP2{}, &gfP2{}, &gfP2{}
	xpoint, dblxpoint, trplxpoint, t1, t2, t3, f := &twistPoint{}, &twistPoint{}, &twistPoint{},
		&twistPoint{}, &twistPoint{}, &twistPoint{}, &twistPoint{}
	for {
		// let's try to construct a point in F(P^2) as 1 + v*i
		x.y = *newGFp(1)
		x.x = *v

		// now we need to check if a is X-coordinate of some point
		// on the curve (if there exists b such that b^2 = a^3 + 3)
		xxx.Square(x)
		xxx.Mul(xxx, x)
		rhs.Add(xxx, twistB)

		y, err := y.Sqrt(rhs)
		if err == nil { // there is a square root for rhs
			point := &twistPoint{
				*x,
				*y,
				gfP2{*newGFp(0), *newGFp(1)},
				gfP2{*newGFp(0), *newGFp(1)},
			}

			// xQ + frob(3*xQ) + frob(frob(xQ)) + frob(frob(frob(Q)))
			// xQ:

			xpoint.Mul(point, u)

			dblxpoint.Double(xpoint)

			trplxpoint.Add(xpoint, dblxpoint)
			trplxpoint.MakeAffine()

			// Frobenius(3*xQ)
			_, err = t1.Frobenius(trplxpoint)
			if err != nil {
				return nil, err
			}

			// Frobenius(Frobenius((xQ))
			xpoint.MakeAffine()
			_, err = t2.Frobenius(xpoint)
			if err != nil {
				return nil, err
			}
			_, err = t2.Frobenius(t2)
			if err != nil {
				return nil, err
			}

			// Frobenius(Frobenius(Frobenius(Q)))
			_, err = t3.Frobenius(point)
			if err != nil {
				return nil, err
			}
			_, err = t3.Frobenius(t3)
			if err != nil {
				return nil, err
			}
			_, err = t3.Frobenius(t3)
			if err != nil {
				return nil, err
			}

			f.Add(xpoint, t1)
			f.Add(f, t2)
			f.Add(f, t3)

			return &G2{f}, nil
		}
		gfpAdd(v, v, newGFp(1))
	}
}

func (e *G2) String() string {
	return "bn256.G2" + e.p.String()
}

// ScalarBaseMult sets e to g*k where g is the generator of the group and then
// returns out.
func (e *G2) ScalarBaseMult(k *big.Int) *G2 {
	if e.p == nil {
		e.p = &twistPoint{}
	}
	e.p.Mul(twistGen, k)
	return e
}

// ScalarMult sets e to a*k and then returns e.
func (e *G2) ScalarMult(a *G2, k *big.Int) *G2 {
	if e.p == nil {
		e.p = &twistPoint{}
	}
	e.p.Mul(a.p, k)
	return e
}

// Add sets e to a+b and then returns e.
func (e *G2) Add(a, b *G2) *G2 {
	if e.p == nil {
		e.p = &twistPoint{}
	}
	e.p.Add(a.p, b.p)
	return e
}

// Neg sets e to -a and then returns e.
func (e *G2) Neg(a *G2) *G2 {
	if e.p == nil {
		e.p = &twistPoint{}
	}
	e.p.Neg(a.p)
	return e
}

// Set sets e to a and then returns e.
func (e *G2) Set(a *G2) *G2 {
	if e.p == nil {
		e.p = &twistPoint{}
	}
	e.p.Set(a.p)
	return e
}

// Marshal converts e into a byte slice.
func (e *G2) Marshal() []byte {
	// Each value is a 256-bit number.
	const numBytes = 256 / 8

	if e.p == nil {
		e.p = &twistPoint{}
	}

	e.p.MakeAffine()
	if e.p.IsInfinity() {
		return make([]byte, 1)
	}

	ret := make([]byte, 1+numBytes*4)
	ret[0] = 0x01
	temp := &gfP{}

	montDecode(temp, &e.p.x.x)
	temp.Marshal(ret[1:])
	montDecode(temp, &e.p.x.y)
	temp.Marshal(ret[1+numBytes:])
	montDecode(temp, &e.p.y.x)
	temp.Marshal(ret[1+2*numBytes:])
	montDecode(temp, &e.p.y.y)
	temp.Marshal(ret[1+3*numBytes:])

	return ret
}

// Unmarshal sets e to the result of converting the output of Marshal back into
// a group element and then returns e.
func (e *G2) Unmarshal(m []byte) ([]byte, error) {
	// Each value is a 256-bit number.
	const numBytes = 256 / 8

	if e.p == nil {
		e.p = &twistPoint{}
	}

	if len(m) > 0 && m[0] == 0x00 {
		e.p.SetInfinity()
		return m[1:], nil
	} else if len(m) > 0 && m[0] != 0x01 {
		return nil, errors.New("bn256: malformed point")
	} else if len(m) < 1+4*numBytes {
		return nil, errors.New("bn256: not enough data")
	}

	e.p.x.x.Unmarshal(m[1:])
	e.p.x.y.Unmarshal(m[1+numBytes:])
	e.p.y.x.Unmarshal(m[1+2*numBytes:])
	e.p.y.y.Unmarshal(m[1+3*numBytes:])
	montEncode(&e.p.x.x, &e.p.x.x)
	montEncode(&e.p.x.y, &e.p.x.y)
	montEncode(&e.p.y.x, &e.p.y.x)
	montEncode(&e.p.y.y, &e.p.y.y)

	if e.p.x.IsZero() && e.p.y.IsZero() {
		// This is the point at infinity.
		e.p.y.SetOne()
		e.p.z.SetZero()
		e.p.t.SetZero()
	} else {
		e.p.z.SetOne()
		e.p.t.SetOne()

		if !e.p.IsOnCurve() {
			return nil, errors.New("bn256: malformed point")
		}
	}

	return m[1+4*numBytes:], nil
}

// GT is an abstract cyclic group. The zero value is suitable for use as the
// output of an operation, but cannot be used as an input.
type GT struct {
	p *gfP12
}

// RandomGT returns x and e(g₁, g₂)ˣ where x is a random, non-zero number read
// from r.
func RandomGT(r io.Reader) (*big.Int, *GT, error) {
	k, err := randomK(r)
	if err != nil {
		return nil, nil, err
	}

	return k, new(GT).ScalarBaseMult(k), nil
}


// returns number in P-representation: a_11*P^11 + ... + a_1*P^1 + a_0 where 0 <= a_i < P
func intToPRepr(n *big.Int) []*big.Int {
	nn := new(big.Int).Set(n)
	pToI := big.NewInt(1)
	mod := new(big.Int).Set(p)
	a := make([]*big.Int, 12)
	for i := 0; i < 12; i++ {
		ai := new(big.Int).Mod(nn, mod)
		nn.Sub(nn, ai)
		ai.Div(ai, pToI)
		a[i] = ai
		if nn.Cmp(big.NewInt(0)) == 0 {
			for {
				i++
				if i == 12 {
					return a
				}
				a[i] = big.NewInt(0)
			}
		}
		pToI.Mul(pToI, p)
		mod.Mul(mod, p)
	}

	return a
}

// converts number in P-representation into *big.Int
func pReprToInt(a []*big.Int) *big.Int {
	pToI := big.NewInt(1)
	n := big.NewInt(0)
	for i := 0; i < 12; i++ {
		t := new(big.Int).Mul(a[i], pToI)
		n.Add(n, t)
		pToI.Mul(pToI, p)
	}

	return n
}

// MapStringToGT maps a string to GT group element. Needed for example when a message to be encrypted
// needs to be mapped into GT group.
func MapStringToGT(msg string) (*GT, error) {
	m := new(big.Int)
	m.SetBytes([]byte(msg))
	bound := new(big.Int).Exp(p, big.NewInt(12), nil)
	if m.Cmp(bound) >= 0 {
		return nil, fmt.Errorf("message is bigger than modulo, use key encapsulation")
	}
	a := intToPRepr(m)
	g := &gfP12{}
	g.x.x.x.SetInt(a[0])
	g.x.x.y.SetInt(a[1])
	g.x.y.x.SetInt(a[2])
	g.x.y.y.SetInt(a[3])
	g.x.z.x.SetInt(a[4])
	g.x.z.y.SetInt(a[5])

	g.y.x.x.SetInt(a[6])
	g.y.x.y.SetInt(a[7])
	g.y.y.x.SetInt(a[8])
	g.y.y.y.SetInt(a[9])
	g.y.z.x.SetInt(a[10])
	g.y.z.y.SetInt(a[11])

	return &GT{g}, nil
}

// MapGTToString maps an element from GT group to a string
func MapGTToString(gt *GT) string {
	a := make([]*big.Int, 12)
	a[0], _ = gt.p.x.x.x.ToInt()
	a[1], _ = gt.p.x.x.y.ToInt()
	a[2], _ = gt.p.x.y.x.ToInt()
	a[3], _ = gt.p.x.y.y.ToInt()
	a[4], _ = gt.p.x.z.x.ToInt()
	a[5], _ = gt.p.x.z.y.ToInt()

	a[6], _ = gt.p.y.x.x.ToInt()
	a[7], _ = gt.p.y.x.y.ToInt()
	a[8], _ = gt.p.y.y.x.ToInt()
	a[9], _ = gt.p.y.y.y.ToInt()
	a[10], _ = gt.p.y.z.x.ToInt()
	a[11], _ = gt.p.y.z.y.ToInt()

	r := pReprToInt(a)
	return string(r.Bytes())
}

// Pair calculates an Optimal Ate pairing.
func Pair(g1 *G1, g2 *G2) *GT {
	return &GT{optimalAte(g2.p, g1.p)}
}

// Miller applies Miller's algorithm, which is a bilinear function from the
// source groups to F_p^12. Miller(g1, g2).Finalize() is equivalent to Pair(g1,
// g2).
func Miller(g1 *G1, g2 *G2) *GT {
	return &GT{miller(g2.p, g1.p)}
}

func (g *GT) String() string {
	return "bn256.GT" + g.p.String()
}

// ScalarBaseMult sets e to g*k where g is the generator of the group and then
// returns out.
func (e *GT) ScalarBaseMult(k *big.Int) *GT {
	if e.p == nil {
		e.p = &gfP12{}
	}
	e.p.Exp(gfP12Gen, k)
	return e
}

// ScalarMult sets e to a*k and then returns e.
func (e *GT) ScalarMult(a *GT, k *big.Int) *GT {
	if e.p == nil {
		e.p = &gfP12{}
	}
	e.p.Exp(a.p, k)
	return e
}

// Add sets e to a+b and then returns e.
func (e *GT) Add(a, b *GT) *GT {
	if e.p == nil {
		e.p = &gfP12{}
	}
	e.p.Mul(a.p, b.p)
	return e
}

// Neg sets e to -a and then returns e.
func (e *GT) Neg(a *GT) *GT {
	if e.p == nil {
		e.p = &gfP12{}
	}
	e.p.Conjugate(a.p)
	return e
}

// Set sets e to a and then returns e.
func (e *GT) Set(a *GT) *GT {
	if e.p == nil {
		e.p = &gfP12{}
	}
	e.p.Set(a.p)
	return e
}

// Finalize is a linear function from F_p^12 to GT.
func (e *GT) Finalize() *GT {
	ret := finalExponentiation(e.p)
	e.p.Set(ret)
	return e
}

// Marshal converts e into a byte slice.
func (e *GT) Marshal() []byte {
	// Each value is a 256-bit number.
	const numBytes = 256 / 8

	if e.p == nil {
		e.p = &gfP12{}
		e.p.SetOne()
	}

	ret := make([]byte, numBytes*12)
	temp := &gfP{}

	montDecode(temp, &e.p.x.x.x)
	temp.Marshal(ret)
	montDecode(temp, &e.p.x.x.y)
	temp.Marshal(ret[numBytes:])
	montDecode(temp, &e.p.x.y.x)
	temp.Marshal(ret[2*numBytes:])
	montDecode(temp, &e.p.x.y.y)
	temp.Marshal(ret[3*numBytes:])
	montDecode(temp, &e.p.x.z.x)
	temp.Marshal(ret[4*numBytes:])
	montDecode(temp, &e.p.x.z.y)
	temp.Marshal(ret[5*numBytes:])
	montDecode(temp, &e.p.y.x.x)
	temp.Marshal(ret[6*numBytes:])
	montDecode(temp, &e.p.y.x.y)
	temp.Marshal(ret[7*numBytes:])
	montDecode(temp, &e.p.y.y.x)
	temp.Marshal(ret[8*numBytes:])
	montDecode(temp, &e.p.y.y.y)
	temp.Marshal(ret[9*numBytes:])
	montDecode(temp, &e.p.y.z.x)
	temp.Marshal(ret[10*numBytes:])
	montDecode(temp, &e.p.y.z.y)
	temp.Marshal(ret[11*numBytes:])

	return ret
}

// Unmarshal sets e to the result of converting the output of Marshal back into
// a group element and then returns e.
func (e *GT) Unmarshal(m []byte) ([]byte, error) {
	// Each value is a 256-bit number.
	const numBytes = 256 / 8

	if len(m) < 12*numBytes {
		return nil, errors.New("bn256: not enough data")
	}

	if e.p == nil {
		e.p = &gfP12{}
	}

	e.p.x.x.x.Unmarshal(m)
	e.p.x.x.y.Unmarshal(m[numBytes:])
	e.p.x.y.x.Unmarshal(m[2*numBytes:])
	e.p.x.y.y.Unmarshal(m[3*numBytes:])
	e.p.x.z.x.Unmarshal(m[4*numBytes:])
	e.p.x.z.y.Unmarshal(m[5*numBytes:])
	e.p.y.x.x.Unmarshal(m[6*numBytes:])
	e.p.y.x.y.Unmarshal(m[7*numBytes:])
	e.p.y.y.x.Unmarshal(m[8*numBytes:])
	e.p.y.y.y.Unmarshal(m[9*numBytes:])
	e.p.y.z.x.Unmarshal(m[10*numBytes:])
	e.p.y.z.y.Unmarshal(m[11*numBytes:])
	montEncode(&e.p.x.x.x, &e.p.x.x.x)
	montEncode(&e.p.x.x.y, &e.p.x.x.y)
	montEncode(&e.p.x.y.x, &e.p.x.y.x)
	montEncode(&e.p.x.y.y, &e.p.x.y.y)
	montEncode(&e.p.x.z.x, &e.p.x.z.x)
	montEncode(&e.p.x.z.y, &e.p.x.z.y)
	montEncode(&e.p.y.x.x, &e.p.y.x.x)
	montEncode(&e.p.y.x.y, &e.p.y.x.y)
	montEncode(&e.p.y.y.x, &e.p.y.y.x)
	montEncode(&e.p.y.y.y, &e.p.y.y.y)
	montEncode(&e.p.y.z.x, &e.p.y.z.x)
	montEncode(&e.p.y.z.y, &e.p.y.z.y)

	return m[12*numBytes:], nil
}
