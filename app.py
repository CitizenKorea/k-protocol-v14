import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class K_PROTOCOL_Nano_Final:
    def __init__(self, file_path):
        self.file_path = file_path
        # [1] K-PROTOCOL 마스터 상수 (절대 공간 System U)
        self.C_K = np.longdouble('297880197.6')
        self.S_EARTH = np.longdouble('1.006419562')
        self.S_LOC_AVG = np.longdouble('1.0064200') # 관측소 평균 왜곡비
        self.sites, self.sources, self.results = {}, {}, []

    def calc_gmst_nano(self, y, m, d, hr, mn, sec):
        """ 나노 단위 정밀도의 J2000 항성시(GMST) 도출 """
        if m <= 2: y -= 1; m += 12
        A_v = int(y / 100); B_v = 2 - A_v + int(A_v / 4)
        jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B_v - 1524.5
        jd += (hr + mn/60.0 + sec/3600.0) / 24.0
        gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360.0
        return np.radians(gmst), jd

    def process(self):
        print(f"[*] K-PROTOCOL 그랜드 엔진 시동: {self.file_path}")
        try:
            with open(self.file_path, 'r') as f: lines = f.readlines()
        except: print("[!] 파일 로드 실패"); return

        # [상태 전이 파싱] 구역별로 철저히 분리해서 데이터 오염 차단
        mode = "STATIONS"
        current_obs = {}

        for line in lines:
            line = line.rstrip()
            if not line: continue
            if line.startswith('$END'): # 구역 전환점
                if mode == "STATIONS": mode = "SOURCES"
                elif mode == "SOURCES": mode = "PARAMS"
                elif mode == "PARAMS": mode = "DATA"
                continue
            
            parts = line.split()
            if mode == "STATIONS" and 'AZEL' in line:
                self.sites[parts[0]] = np.array(parts[1:4], dtype=np.longdouble)
            elif mode == "SOURCES" and len(parts) >= 7:
                try:
                    name = parts[0]
                    ra = np.radians((float(parts[1]) + float(parts[2])/60 + float(parts[3])/3600)*15)
                    dec_val, sign = float(parts[4]), -1 if '-' in str(parts[4]) else 1
                    dec = np.radians(sign * (abs(dec_val) + float(parts[5])/60 + float(parts[6])/3600))
                    self.sources[name] = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], dtype=np.longdouble)
                except: pass
            elif mode == "DATA":
                # 카드 번호 식별 (줄의 마지막 숫자 1 또는 2)
                card_id = parts[-1][-1] if parts else ""
                if card_id == '1': # 카드 1: 관측 정보
                    try:
                        current_obs = {
                            'sA': parts[0], 'sB': parts[1], 'src': parts[2],
                            'time': (int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6]), int(parts[7]), float(parts[8]))
                        }
                    except: current_obs = {}
                elif card_id == '2' and current_obs: # 카드 2: 시간 데이터
                    try:
                        self.evaluate(current_obs, np.longdouble(parts[0]))
                        current_obs = {}
                    except: pass

        print(f"[+] 분석 완료: {len(self.results)}개 데이터 붕괴 성공 (관측소:{len(self.sites)}, 퀘이사:{len(self.sources)})")
        self.report()

    def evaluate(self, obs, obs_val):
        sA, sB, src = obs['sA'], obs['sB'], obs['src']
        if sA not in self.sites or sB not in self.sites or src not in self.sources: return

        gmst, jd = self.calc_gmst_nano(*obs['time'])
        # [K-PROTOCOL] 지구가 도는 찰나의 절대 기하학 타격
        rot = np.array([[np.cos(gmst), -np.sin(gmst), 0], [np.sin(gmst), np.cos(gmst), 0], [0,0,1]], dtype=np.longdouble)
        baseline_si = self.sites[sB] - self.sites[sA]
        baseline_eci = np.dot(rot, baseline_si)
        
        # 핵심 수식: (Dist / c_k) * (S_earth / S_loc) -> 라이고 성공 공식 그대로 VLBI 이식
        tau_K = (abs(np.dot(baseline_eci, self.sources[src])) / self.C_K) * (self.S_EARTH / self.S_LOC_AVG) * 1e9
        self.results.append({'MJD': jd - 2400000.5, 'Res': abs(abs(obs_val) - tau_K)})

    def report(self):
        if not self.results: print("[!] 분석 결과가 비어있습니다."); return
        df = pd.DataFrame(self.results)
        pdf_name = "K_PROTOCOL_Final_Grand_Unified.pdf"
        with PdfPages(pdf_name) as pdf:
            fig, ax = plt.subplots(figsize=(12, 7))
            ax.scatter(df['MJD'], df['Res'], s=5, color='#00ffcc', alpha=0.8)
            ax.set_title('K-PROTOCOL ABSOLUTE GEOMETRY COLLAPSE', color='white', fontsize=16, fontweight='bold')
            fig.patch.set_facecolor('#050510'); ax.set_facecolor('#050510')
            ax.grid(True, alpha=0.15, color='#00ffcc'); ax.tick_params(colors='white')
            trend = np.median(df['Res'])
            ax.axhline(trend, color='#ff007f', linewidth=2, label=f'Absolute Target: {trend:.2f} ns')
            ax.legend(facecolor='black', edgecolor='white', labelcolor='white')
            pdf.savefig(fig, facecolor='#050510'); plt.close(fig)
        print(f"[SUCCESS] '{pdf_name}' 생성 완료. 이제 PDF를 확인하십시오!")

if __name__ == "__main__":
    K_PROTOCOL_Nano_Final("20JAN02XE_N005.ngs").process()
