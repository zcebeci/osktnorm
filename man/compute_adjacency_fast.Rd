\name{compute_adjacency_fast}
\alias{compute_adjacency_fast}
\title{Hızlı Bitişiklik (Adjacency) Matrisi Hesaplama}
\description{
  Ajanların konumlarına ve belirlenen algılama yarıçapına (R) göre sistemin topolojisini belirleyen bitişiklik matrisini oluşturur.
  Tezdeki Denklem 2.2 uyarınca, birbirine R mesafesinden yakın olan ajanlar komşu kabul edilir.
}
\usage{
compute_adjacency_fast(pos, R)
}
\arguments{
  \item{pos}{
    N x n boyutunda sayısal matris. Ajanların öklid uzayındaki konumlarını içerir.
  }
  \item{R}{
    Sayısal değer (skaler). Algılama yarıçapını (Sensing Radius) belirtir.
    İki ajan arasındaki mesafe R'den küçükse aralarında bir bağ olduğu varsayılır.
  }
}
\details{
  Fonksiyon, R'ın yerleşik \code{dist} fonksiyonunu kullanarak tüm ajanlar arasındaki mesafe matrisini çıkarır.
  Ardından vektörel mantıksal karşılaştırma ile şu koşulu kontrol eder:
  \deqn{||x_i - x_j|| < R \quad \textrm{ve} \quad ||x_i - x_j|| > 0}
  Kendisiyle olan bağlantı (mesafe > 0 koşulu ile) hariç tutulur.
  Sonuç, bellek verimliliği için tamsayı (integer) tipine dönüştürülür.
}
\value{
  N x N boyutunda, 1 ve 0'lardan oluşan simetrik bir matris döndürür.
  A(i,j) = 1 ise i ve j ajanları komşudur.
}
\references{
  Chapter 2, Equation 2.2: Neighborhood definition.
  Source Code: Matris Tabanlı Sürü Kontrol Algoritmaları.
}
\author{
  Sürü Kontrol Araştırma Grubu
}
\seealso{
  \code{\link{compute_control_input_fast}}
}
