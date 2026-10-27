"""
NCBI E-utilities için HTTP oturum yöneticisi.

Bu modül, NCBI'a yapılan istekleri rate limiting, otomatik retry ve
connection pooling ile yönetir. ``crn.fetch.crn_fetch_ncbi.Ncbi`` tarafından
kullanılır.
"""

from __future__ import annotations

import time
import threading
from typing import Any, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

class _RateLimiter:
    """İş parçacığı güvenli (thread-safe) hız sınırlayıcı.

    Parameters
    ----------
    min_interval : float
        İki istek arasındaki minimum süre (saniye cinsinden).
    """

    def __init__(self, min_interval: float) -> None:
        self._interval = min_interval
        self._lock     = threading.Lock()
        self._last     = 0.0

    def wait(self) -> None:
        """Gerekirse bir sonraki isteğe kadar bekler."""
        with self._lock:
            elapsed = time.monotonic() - self._last
            gap     = self._interval - elapsed
            if gap > 0:
                time.sleep(gap)
            self._last = time.monotonic()

    @property
    def interval(self) -> float:
        """Mevcut bekleme aralığını döndürür (saniye)."""
        return self._interval

    @interval.setter
    def interval(self, value: float) -> None:
        with self._lock:
            self._interval = value

class NcbiSession:
    """NCBI E-utilities için HTTP oturum yöneticisi.

    API anahtarı verilmezse NCBI'nın ücretsiz limitine (3 istek/sn) uyar;
    API anahtarı verilirse bu limit 10 istek/sn'ye çıkar.

    Parameters
    ----------
    api_key : str, optional
        NCBI API anahtarı. Yoksa anonim limit uygulanır.
    e_mail : str, optional
        NCBI'nın User-Agent için istediği e-posta adresi.
    timeout : int, optional
        Her istek için zaman aşımı süresi (saniye). Varsayılan: 30.
    max_retries : int, optional
        Başarısız isteklerde yeniden deneme sayısı. Varsayılan: 5.
    pool_connections : int, optional
        Bağlantı havuzu sayısı. Varsayılan: 10.
    pool_maxsize : int, optional
        Havuzdaki maksimum bağlantı. Varsayılan: 20.

    Examples
    --------
    Context manager ile kullanım (önerilen):

    >>> with NcbiSession(api_key="...", e_mail="sen@örnek.com") as session:
    ...     response = session.get("https://eutils.ncbi.nlm.nih.gov/...")

    Manuel kullanım:

    >>> session = NcbiSession()
    >>> response = session.get("https://eutils.ncbi.nlm.nih.gov/...")
    >>> session.close()
    """

    _RATE_NO_KEY   = 1 / 3    # 3 istek/sn — API anahtarsız ücretsiz limit
    _RATE_WITH_KEY = 1 / 10   # 10 istek/sn — API anahtarlı limit

    def __init__(
        self,
        api_key:          Optional[str]  = None,
        e_mail:           Optional[str]  = None,
        timeout:          int            = 30,
        max_retries:      int            = 5,
        pool_connections: int            = 10,
        pool_maxsize:     int            = 20,
    ) -> None:
        self._api_key  = api_key
        self._e_mail   = e_mail
        self._timeout  = timeout
        self._closed   = True

        # Rate limiter
        interval       = self._RATE_WITH_KEY if api_key else self._RATE_NO_KEY
        self._limiter  = _RateLimiter(interval)

        # requests.Session kurulumu
        self._session  = requests.Session()
        self._session.headers.update({
            "User-Agent": (
                f"CrnFetchClient/1.0 ({self._e_mail or 'e-posta-yok'})"
            ),
        })
        if api_key:
            self._session.params = {"api_key": api_key}  # type: ignore[assignment]

        retry = Retry(
            total            = max_retries,
            backoff_factor   = 1.2,
            status_forcelist = [429, 500, 502, 503, 504],
            allowed_methods  = ["GET"],
            raise_on_status  = False,
        )
        adapter = HTTPAdapter(
            max_retries      = retry,
            pool_connections = pool_connections,
            pool_maxsize     = pool_maxsize,
        )
        self._session.mount("https://", adapter)
        self._session.mount("http://",  adapter)

        self._closed = False

    def get(self, url: str, **kwargs: Any) -> requests.Response:
        """Rate limiting uygulayarak GET isteği gönderir.

        Parameters
        ----------
        url : str
            İstek yapılacak tam URL.
        **kwargs
            ``requests.Session.get`` metoduna iletilen ek parametreler
            (örn. ``params``, ``headers``).

        Returns
        -------
        requests.Response
            Sunucudan gelen yanıt nesnesi.

        Raises
        ------
        RuntimeError
            Oturum kapalıysa.
        requests.exceptions.RequestException
            Ağ hatalarında veya HTTP hata kodlarında.

        Examples
        --------
        >>> with NcbiSession() as session:
        ...     resp = session.get(
        ...         "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        ...         params={"db": "gene", "term": "BDNF", "retmode": "json"},
        ...     )
        ...     data = resp.json()
        """
        if self._closed:
            raise RuntimeError("Kapalı bir oturum üzerinden istek yapılamaz.")

        self._limiter.wait()
        kwargs.setdefault("timeout", self._timeout)
        response = self._session.get(url, **kwargs)
        response.raise_for_status()
        return response

    def is_open(self) -> bool:
        """Oturumun açık olup olmadığını döndürür.

        Returns
        -------
        bool
            Oturum aktifse ``True``, kapalıysa ``False``.
        """
        return not self._closed

    def close(self) -> None:
        """Oturumu kapatır ve kaynakları serbest bırakır."""
        if not self._closed:
            self._session.close()
            self._closed = True

    def __enter__(self) -> NcbiSession:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def __repr__(self) -> str:
        key_info = f"api_key=***{self._api_key[-4:]}" if self._api_key else "Anonim"
        durum    = "Kapalı" if self._closed else "Açık"
        return (
            f"NcbiSession({key_info}, "
            f"hız={1 / self._limiter.interval:.0f} istek/sn, "
            f"durum={durum})"
        )
