package bn256

import (
	"testing"

	"bytes"
	"crypto/rand"

	"golang.org/x/crypto/bn256"
	"fmt"
	"math"
	"math/big"
)

func TestG1(t *testing.T) {
	k, Ga, err := RandomG1(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(bn256.G1).ScalarBaseMult(k)
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestG1Marshal(t *testing.T) {
	_, Ga, err := RandomG1(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(G1)
	_, err = Gb.Unmarshal(ma)
	if err != nil {
		t.Fatal(err)
	}
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestG2(t *testing.T) {
	k, Ga, err := RandomG2(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(bn256.G2).ScalarBaseMult(k)
	mb := Gb.Marshal()
	mb = append([]byte{0x01}, mb...)

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestG2Marshal(t *testing.T) {
	_, Ga, err := RandomG2(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(G2)
	_, err = Gb.Unmarshal(ma)
	if err != nil {
		t.Fatal(err)
	}
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestGT(t *testing.T) {
	k, Ga, err := RandomGT(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb, ok := new(bn256.GT).Unmarshal((&GT{gfP12Gen}).Marshal())
	if !ok {
		t.Fatal("unmarshal not ok")
	}
	Gb.ScalarMult(Gb, k)
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestGTMarshal(t *testing.T) {
	_, Ga, err := RandomGT(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(GT)
	_, err = Gb.Unmarshal(ma)
	if err != nil {
		t.Fatal(err)
	}
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestBilinearity(t *testing.T) {
	for i := 0; i < 2; i++ {
		a, p1, _ := RandomG1(rand.Reader)
		b, p2, _ := RandomG2(rand.Reader)
		e1 := Pair(p1, p2)

		e2 := Pair(&G1{curveGen}, &G2{twistGen})
		e2.ScalarMult(e2, a)
		e2.ScalarMult(e2, b)

		if *e1.p != *e2.p {
			t.Fatalf("bad pairing result: %s", e1)
		}
	}
}

func TestTripartiteDiffieHellman(t *testing.T) {
	a, _ := rand.Int(rand.Reader, Order)
	b, _ := rand.Int(rand.Reader, Order)
	c, _ := rand.Int(rand.Reader, Order)

	pa, pb, pc := new(G1), new(G1), new(G1)
	qa, qb, qc := new(G2), new(G2), new(G2)

	pa.Unmarshal(new(G1).ScalarBaseMult(a).Marshal())
	qa.Unmarshal(new(G2).ScalarBaseMult(a).Marshal())
	pb.Unmarshal(new(G1).ScalarBaseMult(b).Marshal())
	qb.Unmarshal(new(G2).ScalarBaseMult(b).Marshal())
	pc.Unmarshal(new(G1).ScalarBaseMult(c).Marshal())
	qc.Unmarshal(new(G2).ScalarBaseMult(c).Marshal())

	k1 := Pair(pb, qc)
	k1.ScalarMult(k1, a)
	k1Bytes := k1.Marshal()

	k2 := Pair(pc, qa)
	k2.ScalarMult(k2, b)
	k2Bytes := k2.Marshal()

	k3 := Pair(pa, qb)
	k3.ScalarMult(k3, c)
	k3Bytes := k3.Marshal()

	if !bytes.Equal(k1Bytes, k2Bytes) || !bytes.Equal(k2Bytes, k3Bytes) {
		t.Errorf("keys didn't agree")
	}
}

func TestHashToG1(t *testing.T) {
	s := "hello hello hello hello hello hello hello hello hello hello hello hello hello"
	g1, err := HashG1(s)
	if err != nil {
		t.Errorf("hashing failed: %v", err)
	}
	fmt.Println("hash1", g1)

	x := newGFp(100)
	fmt.Println(x.String())
	fmt.Println(x[0], x[1], x[2], x[3])

	montDecode(x, x)


	xx, bo := new(big.Int).SetString(x.String(), 16)
	fmt.Println("integer", xx, bo)
	z := &gfP{}
	x = newGFp(int64(math.Pow(2, 62)))
	montDecode(z, x)
	fmt.Println(z)

	w := &gfP{}
	gfpMul(w, x, x)
	gfpMul(w, w, w)
	gfpMul(w, w, x)

	montDecode(z, w)

	fmt.Println(z)
	fmt.Println(z[0], z[1], z[2], z[3])
	xx, bo = new(big.Int).SetString(z.String(), 16)
	fmt.Println("integer", xx, bo)

	p, _ := new(big.Int).SetString("65000549695646603732796438742359905742825358107623003571877145026864184071783", 10)
	test := new(big.Int).Exp(big.NewInt(2), big.NewInt(62 * 5), nil)
	test.Mod(test, p)
	fmt.Println("test", test)
	a := test.Text(16)
	fmt.Println(a)

	iii := test.Uint64()
	fmt.Println(iii)
	//aa := uint64(a)
	////strconv

	//i := big.NewInt(5)
	g := intToGfP(test)
	//gh := &gfP{}
	//montDecode(gh, g)
	fmt.Println(g)
	fmt.Println(g[0], g[1], g[2], g[3])
	montDecode(z, g)

	fmt.Println(z)
	fmt.Println(z[0], z[1], z[2], z[3])

	test2, err := gfPToInt(g)
	fmt.Println(test2, err)

}

func TestHashToG2(t *testing.T) {
	s := "hello hello hello hello hello hello hello hello hello hello hello hello hello"
	g2, err := HashG2(s)
	if err != nil {
		t.Errorf("hashing failed: %v", err)
	}
	fmt.Println(g2)
}


func BenchmarkG1(b *testing.B) {
	x, _ := rand.Int(rand.Reader, Order)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		new(G1).ScalarBaseMult(x)
	}
}

func BenchmarkG2(b *testing.B) {
	x, _ := rand.Int(rand.Reader, Order)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		new(G2).ScalarBaseMult(x)
	}
}

func BenchmarkGT(b *testing.B) {
	x, _ := rand.Int(rand.Reader, Order)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		new(GT).ScalarBaseMult(x)
	}
}

func BenchmarkPairing(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Pair(&G1{curveGen}, &G2{twistGen})
	}
}
